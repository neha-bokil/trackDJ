#' @importFrom methods is
#' @importFrom stats start end
#' @importFrom utils read.delim read.table data
#' @importFrom rlang .data
#' @importFrom data.table :=
#' @importFrom ggplot2 facet_grid
NULL

#' plot gene models
#' @param genomicLoc character vector. A gene name (ex:"ZFX") OR a genomic locus of the form c(chromosome, start coordinate, end coordinate) (ex: c("chrX", 24040226, 24232664))
#' @param mart a mart obtained by utilizing the useMart function in biomaRt
#' @param ensembl_set if 'mart' is not provided, which Ensembl annotation you want to pull from. Default is 'hsapiens_gene_ensembl'. Most ensembl sets follow this format; for example if you're using mouse it would be 'mmusculus_gene_ensembl'
#' @param gene_symbol nomenclature system you are using for genes. Default is 'hgnc_symbol'but switches to 'external_gene_name' for any ensembl_set that is not 'hsapiens_gene_ensembl'. If you are using other nomenclature systems (ex: 'mgi_symbol' for mouse), you must specify. If using a custom annotation file, you likely need to use 'gene_id'.
#' @param includeNPC boolean, whether to include genomic features that are not protein coding. Default is FALSE.
#' @param custom_anno file path to a .gtf/.gff if using a custom annotation rather than ensembl
#' @param includeTranscripts boolean. Whether to plot transcripts. Default is FALSE.
#' @param includeTxtNames boolean. Whether to label transcript models. If set to TRUE, requires includeGenome and include Transcripts to be set to TRUE. Default is FALSE.
#' @param transcript_list list, used to specify the names of the transcripts you want to include (if you only want a specific set of transcripts). Requires includeTranscripts to be set to TRUE. overrides any specified transcript filters
#' @param supportedTranscriptsOnly boolean, whether to only include well-supported transcripts as defined by the 'transcript_filters' parameter. Default is TRUE.
#' @param transcript_filters character vector of which transcript filters to clear when supportedTranscriptsOnly is set to TRUE. Common choices from ensembl are 'transcript_gencode_basic', 'transcript_appris', and 'transcript_tsl'. the 'appris_options' and 'tsl_options' parameters allow you to specify which APPRIS classification or transcript support levels, respectively, that you want to keep. Any other filter will simply keep transcripts that have a non-empty field for that filter. Default is 'transcript_gencode_basic'
#' @param appris_options character vector of APPRIS filters to implement if 'transcript_appris' is provided in 'transcript_filters'. can specify both principal (P1-P5) and alternative transcripts (A1-A5). ex: c("P1","P2","A1") will give you principal 1, principal2, and alternative 1 transcripts. If 'transcript_appris' is provided as 'transcript_filter' and 'appris_options' is not specified, all principal transcripts will be kept
#' @param tsl_options character vector of transcript support levels you want to keep if 'transcript_tsl' is provided in 'transcript_filters'
#' @param tag_options character vector of which tags to keep if using a custom annotation and setting supportedTranscriptsOnly to TRUE.
#' @param canonicalTranscriptOnly boolean, whether to only show the canoncial transcript. Default is FALSE.
#' @param upDown a character vector of length 2. if genomicLoc is a gene, how many base pairs upstream (first value) and downstream (second value) of the gene you want in your figure. Default is c(2000,2000)
#' @param fontSize a numeric for desired font size. Default is 9.
#' @return list of length 2: figure (ggplot of gene models); genePlot (ggplot of gene models)

plot_gene<-function(genomicLoc, includeTranscripts=FALSE,includeTxtNames = TRUE, transcript_list=NULL,
                    supportedTranscriptsOnly=TRUE, transcript_filters=c("transcript_gencode_basic"),
                    appris_options=NULL,tsl_options=NULL,tag_options=NULL,
                    canonicalTranscriptOnly=FALSE, mart=NULL,
                    ensembl_set="hsapiens_gene_ensembl",
                    gene_symbol="hgnc_symbol", includeNPC=FALSE, custom_anno=NULL,upDown=c(2000,2000), fontSize=9) {
  fontSize_text<-fontSize/ggplot2::.pt
  #genomicLoc should either be in the form of a gene name (string) or as a vector with three entries (chromosome, start coordinate, end coordinate)
  if (length(genomicLoc) !=1 & length(genomicLoc)!=3) {
    stop("error: Please provide either a gene name or set of genomic coordinates")
  }

  #check that only one of mart or custom anno is used
  if(is.null(mart)==FALSE & is.null(custom_anno)==FALSE) {
    stop("error: Please choose one of mart or custom_anno")
  }

  #when no custom annotation is provided, need to pull from ensembl mart
  if (is.null(custom_anno)==TRUE) {
    if (is.null(mart)==TRUE) {
      #no mart provided - need to set it up
      if (ensembl_set != "hsapiens_gene_ensembl") {
        #cannot use "hgnc_symbol" when not working in human - switch to "external_gene_name"
        if (gene_symbol == "hgnc_symbol") {
          gene_symbol<-"external_gene_name"
          message("cannot use hgnc_symbol for non-human cases. switching to external_gene_name")
        }
      }
      #set up mart
      mart<-biomaRt::useMart("ensembl", dataset=ensembl_set)
    }
    else {
      ensembl_set<-mart@dataset
    }
    #make sure gene_symbol is in mart, otherwise switch to external_gene_name
    if (gene_symbol %in% biomaRt::listAttributes(mart)$name == FALSE) {
      if ("external_gene_name" %in% biomaRt::listAttributes(mart)$name==TRUE) {
        gene_symbol<-"external_gene_name"
        message("provided gene_symbol not found in mart. switching to external_gene_name")
      }
      else {
        stop("error: please provide valid gene_symbol present in mart")
      }
    }

    #if genomicLoc is a gene:
    if (length(genomicLoc)==1) {
      geneName<-genomicLoc
      #get coordinate information of gene

      gene_info<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position","strand"),filters = gene_symbol,values = geneName,mart = mart)}, error=function(e) {
        message("getBM failed: ", e$message)
        NULL
      })

      if (is.null(gene_info)) {
        message("trying useast mirror")
        eastMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "useast")
        gene_info<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position","strand"),filters = gene_symbol,values = geneName,mart = eastMart)}, error=function(e) {
          message("getBM with useast mirror failed: ", e$message)
          NULL
        })
        if (is.null(gene_info)) {
          message("trying www mirror")
          westMart<-biomaRt::useEnsembl(biomart="genes", dataset=ensembl_set, mirror="www")
          gene_info<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position","strand"),filters = gene_symbol,values = geneName,mart = westMart)}, error=function(e) {
            message("getBM with www mirror failed: ", e$message)
            NULL
          })

          if (is.null(gene_info)) {
            message("trying asia mirror")
            asiaMart<-biomaRt::useEnsembl(biomart="genes", dataset=ensembl_set, mirror="asia")
            gene_info<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position","strand"),filters = gene_symbol,values = geneName,mart = asiaMart)}, error=function(e) {
              message("getBM with asia mirror failed: ", e$message)
              NULL
            })
            if (is.null(gene_info)) {
              stop("error: all biomaRt::getBM() attempts failed")
            }

          }

          }
      }

      if (base::nrow(gene_info)>1) {
        message("warning: more than one location detected for the provided gene name")
        if (("X" %in% gene_info$chromosome_name | "chrX" %in% gene_info$chromosome_name) & ("Y" %in% gene_info$chromosome_name | "chrY" %in% gene_info$chromosome_name)) {
          message("Gene is present on both X and Y chromosomes - only chrX coordinates will be used. please use coordinates if looking at chrY")
          gene_info<-gene_info[gene_info$chromosome_name=="X" | gene_info$chromosome_name=="chrX",]
          if (base::nrow(gene_info)>1) {
            message("warning: more than one location detected on chrX - using the first. if you meant a different location, please use coordinates instead")
          }
        }
        if (("Z" %in% gene_info$chromosome_name | "chrZ" %in% gene_info$chromosome_name) & ("W" %in% gene_info$chromosome_name | "chrW" %in% gene_info$chromosome_name)) {
          message("Gene is present on both Z and W chromosomes - only chrZ coordinates will be used. please use coordinates if looking at chrW")
          gene_info<-gene_info[gene_info$chromosome_name=="Z" | gene_info$chromosome_name=="chrZ",]
          if(base::nrow(gene_info)>1) {
            message("warning: more than one location detected on chrZ - using the first. if you meant a different location, please use coordinates instead")

          }
        }

        gene_info<-gene_info[1,]

      }

      #up stream and downstream directions will depend on strand
      if (is.na(gene_info$strand[1])==FALSE & (gene_info$strand[1]=="+" | gene_info$strand[1]==1)) {
        #gene is on + strand
        start <- gene_info$start_position[1] - upDown[1]
        end <- gene_info$end_position[1] + upDown[2]

      }
      else if (is.na(gene_info$strand[1])==FALSE & (gene_info$strand[1]=="-" | gene_info$strand[1]==(-1))){
        #gene is on - strand
        start <- gene_info$end_position[1] + upDown[1]
        end <- gene_info$start_position[1] - upDown[2]
      }
      else {
        message("strand info not provided. assuming + strand")
        start <- gene_info$start_position[1] - upDown[1]
        end <- gene_info$end_position[1] + upDown[2]
      }

      chr <- gene_info$chromosome_name[1]
      gene_info$chromosome_name<-paste0("chr", gene_info$chromosome_name)

      graphTitle<-paste0(genomicLoc, "\n", gene_info$chromosome_name, ":",min(as.numeric(start), as.numeric(end)), "-", max(as.numeric(start),as.numeric(end))) #used for the title of the plot
    }


    else {
      #if coordinates are provided rather than a gene name
      upDown<-c(0,0)
      geneName<-"none"

      #provided locus may or may not include 'chr' when indicating chromosome - need to account for both
      if (grepl("chr", genomicLoc[1])==TRUE) {
        graphTitle<-paste0(genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
        gene_info<-data.frame(chromosome_name=genomicLoc[1], start_position=genomicLoc[2], end_position=genomicLoc[3])
        chr <- substr(genomicLoc[1], 4, nchar(genomicLoc[1]))
      }
      else {
        graphTitle<-paste0("chr",genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
        gene_info<-data.frame(chromosome_name=paste0("chr",genomicLoc[1]), start_position=genomicLoc[2], end_position=genomicLoc[3])
        chr<-genomicLoc[1]
      }

      gene_info$start_position<-as.numeric(gene_info$start_position)
      gene_info$end_position<-as.numeric(gene_info$end_position)

      start <- as.numeric(gene_info$start_position[1])
      end <- as.numeric(gene_info$end_position[1])

    }

    #grab gene annotations. will need "gene_biotype" to filter out non-proteinCoding genes, so ensure it is an available attribute
    if ("gene_biotype" %in% biomaRt::listAttributes(mart)$name==FALSE) {
      typeFilter<-NULL
      if (includeNPC==FALSE) {
        warning("'gene_biotype' unavailable in mart. unable to filter out non-protein-coding genes")
      }
    }
    else {
      typeFilter<-"gene_biotype"
    }

    martUsed<-"og"
    gene_annot<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", gene_symbol, "strand", typeFilter),
                                         filters = c("chromosome_name", "start", "end"),
                                         values = list(chr, min(as.numeric(start),as.numeric(end)), max(as.numeric(start),as.numeric(end))),
                                         mart = mart)}, error=function(e) {
      message("getBM failed: ", e$message)
      NULL
    })

    if (is.null(gene_annot)) {
      message("trying useast mirror")
      martUsed<-"east"
      if (exists("eastMart")==FALSE) {
        eastMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "useast")
      }

      gene_annot<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", gene_symbol, "strand", typeFilter),
                                           filters = c("chromosome_name", "start", "end"),
                                           values = list(chr, min(as.numeric(start),as.numeric(end)), max(as.numeric(start),as.numeric(end))),
                                           mart = eastMart)}, error=function(e) {
                                             message("getBM with useast mirror failed: ", e$message)
                                             NULL
                                           })
      if (is.null(gene_annot)) {
        message("trying www mirror")
        martUsed<-"west"
        if (exists("westMart")==FALSE) {
          westMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "www")
        }

        gene_annot<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", gene_symbol, "strand", typeFilter),
                                             filters = c("chromosome_name", "start", "end"),
                                             values = list(chr, min(as.numeric(start),as.numeric(end)), max(as.numeric(start),as.numeric(end))),
                                             mart = westMart)}, error=function(e) {
                                               message("getBM with www mirror failed: ", e$message)
                                               NULL
                                             })

        if (is.null(gene_annot)) {
          message("trying asia mirror")
          martUsed<-"asia"
          if (exists("asiaMart")==FALSE) {
            asiaMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "asiaMart")
          }

          gene_annot<-tryCatch({biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", gene_symbol, "strand", typeFilter),
                                               filters = c("chromosome_name", "start", "end"),
                                               values = list(chr, min(as.numeric(start),as.numeric(end)), max(as.numeric(start),as.numeric(end))),
                                               mart = asiaMart)}, error=function(e) {
                                                 message("getBM with asia mirror failed: ", e$message)
                                                 NULL
                                               })

          if (is.null(gene_annot)) {
            martUsed<-"none"
            stop("error: all biomaRt::getBM() attempts failed")
          }
        }

      }
    }


    #filter out non-protein-coding genes, if applicable
    if (includeNPC==FALSE) {
      #need gene_biotype information for filtering
      if ("gene_biotype" %in% colnames(gene_annot)) {
        if (length(genomicLoc)==1) {
          #do not want to remove given gene even if it is not protein coding
          toExclude<-gene_annot[gene_annot$gene_biotype != "protein_coding" & gene_annot[[gene_symbol]] != genomicLoc,][[gene_symbol]]
        }
        else {
          toExclude<-gene_annot[gene_annot$gene_biotype != "protein_coding",][[gene_symbol]]
        }
        gene_annot<-gene_annot[gene_annot[[gene_symbol]] %in% toExclude ==FALSE,]
      }
    }

    if (base::nrow(gene_annot)!=0) {
      names(gene_annot)[names(gene_annot)=="start_position"]<-"start"
      names(gene_annot)[names(gene_annot)=="end_position"]<-"end"
      gene_annot$type<-"Gene"
      gene_annot$label<-gene_annot[[gene_symbol]]

      #order the genes - if genomicLoc was a gene name, it should be first
      if (nrow(gene_annot[gene_annot$gene ==geneName,])!=0) {
        if (nrow(gene_annot[gene_annot$gene != geneName,])!=0) {
          otherGenes<-gene_annot[gene_annot$gene != geneName,]
          otherGenes<-otherGenes[order(otherGenes$start),]
          gene_annot<-rbind(gene_annot[gene_annot$gene==geneName,], otherGenes)
        }
      }
      else {
        gene_annot<-gene_annot[order(gene_annot$start),]
      }

      gene_annot<-assign_y_levels_genes(gene_annot)


      gene_annot$trackColor<-"black"
      arrow_properties_gene <- grid::arrow(type = "closed", angle = 20, length = grid::unit(0.3, "cm"))
      gene_annot$rowNum<-1:base::nrow(gene_annot)

      gene_annot$start_real<-gene_annot$start
      gene_annot$end_real<-gene_annot$end

      for(geneNum in gene_annot$rowNum) {
        #plotting genes with start/end outside of range may result in error. adjust gene coordinates to fall in range, if necessary:
        startCoord<-gene_annot[gene_annot$rowNum==geneNum,]$start_real
        endCoord<-gene_annot[gene_annot$rowNum==geneNum,]$end_real
        strand<-gene_annot[gene_annot$rowNum==geneNum,]$strand
        gene_annot[gene_annot$rowNum==geneNum,]$start<-max(startCoord,min(as.numeric(start),as.numeric(end)))
        gene_annot[gene_annot$rowNum==geneNum,]$end<-min(endCoord,max(as.numeric(start),as.numeric(end)))

        #arrow direction depends on whether a gene is on + strand or - strand. start and end should be reversed if gene is on - strand
        gene_annot2<-gene_annot
        if (is.na(strand)==FALSE & (strand==(-1)| strand=="-")) {
          gene_annot[gene_annot$rowNum==geneNum,]$start<-gene_annot2[gene_annot2$rowNum==geneNum,]$end
          gene_annot[gene_annot$rowNum==geneNum,]$end<-gene_annot2[gene_annot2$rowNum==geneNum,]$start

        }
      }

      gene_annot$strand<-as.character(gene_annot$strand)
    }





  }

  else { #what to do with custom annotation
    if (is.character(custom_anno) == TRUE) {
      #check that the annotation file exists and import it if it does
      if (file.exists(custom_anno) == TRUE) {
        message("importing custom annotation...")
        annotation<-rtracklayer::import(custom_anno)
        message("done!")
      }
      else {
        stop("error: custom annotation file does not exist")
      }

    }
    else {
      #if custom_anno is not a file, it should be a GRanges object
      if (is(custom_anno,"GRanges") == FALSE) {
        stop("error: custom annotation must be a gtf, gff, or GRanges object")
      }
      else {
        annotation<-custom_anno
      }
    }

    #check that gene_symbol is in the annotation
    if (gene_symbol %in% colnames(as.data.frame(annotation))==FALSE) {
      if ("gene_name" %in% colnames(as.data.frame(annotation))==TRUE) {
        gene_symbol<-"gene_name"
        message("provided gene_symbol not found in provided annotation. switching to gene_name")
      }
      else {
        if ("gene_id" %in% colnames(as.data.frame(annotation))==TRUE) {
          gene_symbol <- "gene_id"
          message("provided gene_symbol not found in provided annotation. switching to gene_id")
        }
        else {
          stop("error: provided gene_symbol not found in provided annotation")
        }
      }
    }

    #keep only genes
    genes<-as.data.frame(annotation[annotation$type=="gene"])
    #if genomicLoc is a gene:
    if (length(genomicLoc)==1) {
      geneName<-genomicLoc
      gene_info <- genes[genes[[gene_symbol]] == genomicLoc,c("seqnames","start","end","strand")]
      if (base::nrow(gene_info)>1) {
        message("warning: more than one location detected for the provided gene name")
        if (("X" %in% gene_info$seqnames | "chrX" %in% gene_info$seqnames) & ("Y" %in% gene_info$seqnames | "chrY" %in% gene_info$seqnames)) {
          message("Gene is present on both X and Y chromosomes - only chrX coordinates will be used. please use coordinates if looking at chrY")
          gene_info<-gene_info[gene_info$seqnames=="X" | gene_info$seqnames=="chrX",]
          if (base::nrow(gene_info)>1) {
            message("warning: more than one location detected on chrX - using the first. if you meant a different location, please use coordinates instead")
          }
        }
        if (("Z" %in% gene_info$seqnames | "chrZ" %in% gene_info$seqnames) & ("W" %in% gene_info$seqnames | "chrW" %in% gene_info$seqnames)) {
          message("Gene is present on both Z and W chromosomes - only chrZ coordinates will be used. please use coordinates if looking at chrW")
          gene_info<-gene_info[gene_info$seqnames=="Z" | gene_info$seqnames=="chrZ",]
          if(base::nrow(gene_info)>1) {
            message("warning: more than one location detected on chrZ - using the first. if you meant a different location, please use coordinates instead")

          }
        }

        gene_info<-gene_info[1,]
      }
      chr <- gene_info$seqnames

      #need to know strand to determine upstream and downstream directions
      if (is.na(gene_info$strand)==FALSE & (gene_info$strand=="+"|gene_info$strand==1)) {
        start <- gene_info$start[1] - upDown[1]
        end <- gene_info$end[1] + upDown[2]

      }
      else if (is.na(gene_info$strand)==FALSE & (gene_info$strand=="-" | gene_info$strand==(-1))) {
        start <- gene_info$start[1] + upDown[1]
        end <- gene_info$end[1] - upDown[2]
      }
      else {
        message("strand info not provided. assuming + strand")
        gene_info$strand<-"+"
        start <- gene_info$start[1] - upDown[1]
        end <- gene_info$end[1] + upDown[2]
      }

      if (grepl("chr", as.character(chr))==FALSE) {
        if (suppressWarnings(is.na(as.numeric(as.character(chr))))==FALSE | as.character(chr) %in% c("X","Y","Z","W")) {
          gene_info$chromosome_name<-paste0("chr", gene_info$chromosome_name)
          chr<-paste0("chr",as.character(chr))
        }
      }


      graphTitle<-paste0(genomicLoc, "\n", gene_info$seqnames, ":",min(as.numeric(start),as.numeric(end)), "-", max(as.numeric(start),as.numeric(end))) #used for the title of the plot

    }

    else {
      #if coordinates provided
      geneName<-"none"
      upDown<-c(0,0)
      #provided locus may or may not include 'chr' when indicating chromosome - need to account for both


      if (grepl("chr", genomicLoc[1])==TRUE) {
        graphTitle<-paste0(genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
        gene_info<-data.frame(chromosome_name=genomicLoc[1], start_position=genomicLoc[2], end_position=genomicLoc[3])
        chr <- genomicLoc[1]
      }
      else {
        if (suppressWarnings(is.na(as.numeric(as.character(genomicLoc[1]))))==FALSE | as.character(genomicLoc[1]) %in% c("X","Y","Z","W")) {
          graphTitle<-paste0("chr",genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
          gene_info<-data.frame(chromosome_name=paste0("chr",genomicLoc[1]), start_position=genomicLoc[2], end_position=genomicLoc[3])
          chr<-paste0("chr",genomicLoc[1])
        }
        else {
          graphTitle<-paste0(genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
          gene_info<-data.frame(chromosome_name=genomicLoc[1], start_position=genomicLoc[2], end_position=genomicLoc[3])
          chr <- genomicLoc[1]
        }

      }
      gene_info$start<-as.numeric(gene_info$start)
      gene_info$end<-as.numeric(gene_info$end)

      start <- gene_info$start[1]
      end <- gene_info$end[1]
    }


    #find all genes associated with locus

    genes_gr<-annotation[annotation$type=="gene"]

    region_gr <- GenomicRanges::GRanges(seqnames = chr,ranges = IRanges::IRanges(start = min(as.numeric(start),as.numeric(end)), end = max(as.numeric(start),as.numeric(end))))

    overlap_gr <- GenomicRanges::findOverlaps(genes_gr, region_gr)

    genesInRegion <- genes_gr[S4Vectors::queryHits(overlap_gr)]


    gene_annot <- data.frame( chr = as.character(GenomeInfoDb::seqnames(genesInRegion)),
                              gene = as.data.frame(genesInRegion)[[gene_symbol]],
                              start = IRanges::start(genesInRegion),
                              end = IRanges::end(genesInRegion),
                              strand = GenomicRanges::strand(genesInRegion))






    #exclude non-protein-coding genes, if applicable
    if (includeNPC==FALSE) {
      if ("gene_type" %in% colnames(as.data.frame(genesInRegion)) ==FALSE &"gene_biotype" %in% colnames(as.data.frame(genesInRegion))==FALSE) {
        warning("'gene_type' and 'gene_biotype' unavailable in provided annotation. unable to filter out non-protein-coding genes")
      }

      else {
        if ("gene_type" %in% colnames(as.data.frame(genesInRegion) )==FALSE & "gene_biotype" %in% colnames(as.data.frame(genesInRegion))==TRUE) {
          gene_annot$gene_biotype<-as.data.frame(genesInRegion)$gene_biotype
        }
        else {
          gene_annot$gene_biotype<-as.data.frame(genesInRegion)$gene_type
        }

        if (length(genomicLoc)==1) { #exclude all non-protein-coding genes, but if genomicLoc is a non-protein-coding gene, keep it
          toExclude<-gene_annot[gene_annot$gene_biotype != "protein_coding" & gene_annot$gene != genomicLoc,]$gene
        }
        else { #if coordinates are provided, exclude all non-protein-coding genes
          toExclude<-gene_annot[gene_annot$gene_biotype != "protein_coding",]$gene
        }
        gene_annot<-gene_annot[gene_annot$gene %in% toExclude ==FALSE,]

      }


    }

    if (base::nrow(gene_annot)!=0) {
      gene_annot$type<-"Gene"
      gene_annot$label<-gene_annot$gene


      ######
      if (nrow(gene_annot[gene_annot$gene ==geneName,])!=0) {
        if (nrow(gene_annot[gene_annot$gene != geneName,])!=0) {
          otherGenes<-gene_annot[gene_annot$gene != geneName,]
          otherGenes<-otherGenes[order(otherGenes$start),]
          gene_annot<-rbind(gene_annot[gene_annot$gene==geneName,], otherGenes)
        }
      }
      else {
        gene_annot<-gene_annot[order(gene_annot$start),]
      }

      gene_annot<-assign_y_levels_genes(gene_annot)

      gene_annot$trackColor<-"black"
      arrow_properties_gene <- grid::arrow(type = "closed", angle = 20, length = grid::unit(0.3, "cm"))
      gene_annot$rowNum<-1:base::nrow(gene_annot)


      gene_annot$start_real<-gene_annot$start
      gene_annot$end_real<-gene_annot$end
      for(geneNum in gene_annot$rowNum) {
        #plotting genes with start/end outside of range may result in error. adjust gene coordinates to fall in range, if necessary:
        startCoord<-gene_annot[gene_annot$rowNum==geneNum,]$start_real
        endCoord<-gene_annot[gene_annot$rowNum==geneNum,]$end_real
        strand<-gene_annot[gene_annot$rowNum==geneNum,]$strand
        gene_annot[gene_annot$rowNum==geneNum,]$start<-max(startCoord,min(as.numeric(start),as.numeric(end)))
        gene_annot[gene_annot$rowNum==geneNum,]$end<-min(endCoord,max(as.numeric(start),as.numeric(end)))

        #take into account strand information for arrow direction
        gene_annot2<-gene_annot
        if (is.na(strand)==FALSE & (strand==(-1)| strand=="-")) {
          gene_annot[gene_annot$rowNum==geneNum,]$start<-gene_annot2[gene_annot2$rowNum==geneNum,]$end
          gene_annot[gene_annot$rowNum==geneNum,]$end<-gene_annot2[gene_annot2$rowNum==geneNum,]$start

        }
      }


      gene_annot$strand<-as.character(gene_annot$strand)

      names(gene_annot)[names(gene_annot)=="gene"]<-gene_symbol

    }

  }

  #if we aren't plotting transcripts, we can make the plot now
  if (includeTranscripts == FALSE) {
    if (base::nrow(gene_annot)!=0) {

      gene_annot$start<-gene_annot$start/1000000
      gene_annot$end<-gene_annot$end/1000000

      gene_annot$midX<-(gene_annot$start+gene_annot$end)/2
      gene_annot$trackType<-"gene"
      gene_annot$plot_start<-min(start,end)
      gene_annot$plot_end<-max(start,end)

      p1<-ggplot2::ggplot(gene_annot, ggplot2::aes(x = .data$start, y = .data$graph_location)) +
        ggplot2::geom_segment(ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                     data=gene_annot, linewidth=2.5, arrow=arrow_properties_gene, lineend = "butt",
                     linejoin = "mitre")+
        ggplot2::geom_label(color="black", fill="white",ggplot2::aes(x=.data$midX, y=.data$graph_location, label=.data$label),
                   fontface = "italic",size = fontSize_text)+
        ggplot2::xlim(min(as.numeric(start),as.numeric(end))/1000000, max(as.numeric(start),as.numeric(end))/1000000)+
        ggplot2::ylim(min(gene_annot$graph_location-2),2)+ggplot2::theme_classic()+ ggplot2::ggtitle(graphTitle)+
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              axis.line.y = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5))+
        ggplot2::xlab("Location (Mb)")
    }
    else { #no genes in region - blank plot
      gene_annot<-data.frame(chromosome_name=chr, start=as.numeric(start)/1000000, end=as.numeric(end)/1000000, gene_id=NA, strand=NA,
                             type=NA, label=NA, trackColor=NA, rowNum=NA, start_real=start, end_real=end, graph_location=0, midX=((as.numeric(start)+as.numeric(end))/2)/1000000, trackType="gene", plot_start=start, plot_end=end )
      p1<-ggplot2::ggplot(gene_annot)+ggplot2::geom_blank() + ggplot2::xlim(min(as.numeric(start),as.numeric(end))/1000000, max(as.numeric(start),as.numeric(end))/1000000) +
        ggplot2::ylim(1,2) +ggplot2::theme_classic()+ ggplot2::ggtitle(graphTitle)+
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              axis.line.y = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5))+
        ggplot2::xlab("Location (Mb)")

    }

  }

  #getting transcript info, if required
  else { #what to do if includeTranscripts=TRUE

    #get transcript annotations from ensembl when not using a custom annotation
    if (is.null(custom_anno)==TRUE) {
      gene_annot_state<-FALSE
      if (base::nrow(gene_annot) !=0) {
        gene_annot$displayName<-gene_annot[[gene_symbol]]
        gene_annot_state<-TRUE

      }


      #only use available transcript filters
      transcript_filters_present<-transcript_filters[transcript_filters %in% biomaRt::listAttributes(mart)$name==TRUE]
      if ((length(transcript_filters_present) < length(transcript_filters)) & is.null(transcript_list)==FALSE) {
        message(paste(c("transcript filters not found:", transcript_filters[transcript_filters%in% biomaRt::listAttributes(mart)$name==FALSE],"- skipping"), collapse=" "))
      }
      #from mart, "transcript support level" is usually "transcript_tsl"
      if (("transcript_support_level" %in% transcript_filters==TRUE) & ("transcript_support_level" %in% transcript_filters_present ==FALSE ) & ("transcript_tsl" %in% transcript_filters_present == FALSE) & ("transcript_tsl" %in% biomaRt::listAttributes(mart)$name == TRUE)) {
        transcript_filters_present<-c(transcript_filters_present,"transcript_tsl")
      }


      #get transcript annotations


      transcript_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "strand", gene_symbol,"external_transcript_name","transcript_is_canonical",transcript_filters_present),
                                                 filters = c("chromosome_name", "start", "end"),
                                                 values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                                 mart = mart)}, error=function(e) {
        message("getBM failed: ", e$message)
        NULL
      })

      if (is.null(transcript_annot)) {
        message("trying useast mirror")
        if (exists("eastMart")==FALSE) {
          eastMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "useast")
        }

        transcript_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "strand", gene_symbol,"external_transcript_name","transcript_is_canonical",transcript_filters_present),
                                                   filters = c("chromosome_name", "start", "end"),
                                                   values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                                   mart = eastMart)}, error=function(e) {
                                                     message("getBM failed: ", e$message)
                                                     NULL
                                                   })

        if (is.null(transcript_annot)) {
          message("trying www mirror")
          if (exists("westMart")==FALSE) {
            westMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "www")
          }
          transcript_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "strand", gene_symbol,"external_transcript_name","transcript_is_canonical",transcript_filters_present),
                                                     filters = c("chromosome_name", "start", "end"),
                                                     values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                                     mart = westMart)}, error=function(e) {
                                                       message("getBM failed: ", e$message)
                                                       })
          if (is.null(transcript_annot)) {
            message("trying asia mirror")
            if (exists("asiaMart")==FALSE) {
              asiaMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "asia")
            }

            transcript_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "strand", gene_symbol,"external_transcript_name","transcript_is_canonical",transcript_filters_present),
                                                       filters = c("chromosome_name", "start", "end"),
                                                       values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                                       mart = asiaMart)}, error=function(e) {
                                                         message("getBM failed: ", e$message)
                                                       })

            if (is.null(transcript_annot)) {
              stop("failed to get transcript annotation")
            }
          }
        }
      }

      transcript_annot<-data.frame(transcript_annot, stringsAsFactors = FALSE, check.names=FALSE)

      if ("transcript_support_level" %in% colnames(transcript_annot)==TRUE) {
        names(transcript_annot)[names(transcript_annot)=="transcript_support_level"]<-"transcript_tsl"

      }

      #filtering transcripts
      transcript_annot<-transcript_annot[transcript_annot[[gene_symbol]] %in% gene_annot[[gene_symbol]] | transcript_annot$ensembl_transcript_id %in% transcript_list ,] #only use transcripts associated with the genes in gene_annot or provided in transcript_list

      if (is.null(transcript_list)==FALSE) { #if a list of transcripts is provided, keep only those and ignore any specified filters
        transcript_annot<-transcript_annot[transcript_annot$ensembl_transcript_id %in% transcript_list,]
        if (base::nrow(transcript_annot)==0) {
          message("no transcripts found from provided list. will not display any transcript info")
          transcript_annot<-NULL
        }
        else if (base::nrow(transcript_annot)<length(transcript_list)) {
          message(paste(c("transcripts not found:",transcript_list[transcript_list %in% transcript_annot$ensembl_transcript_id == FALSE], "- ignoring" ), collapse=" "))
        }
      }
      else if (supportedTranscriptsOnly==TRUE & canonicalTranscriptOnly==FALSE) { #if no list of transcripts provided and need to use provided filters to get only supported transcripts
        to_keep_all<-c()
        for (filterName in transcript_filters_present) {
          #keep transcripts that pass each filter
          if (filterName != "transcript_appris" & filterName !="transcript_tsl") {
            #for these filters, just want each transcript with a non-empty field for that filter
            to_keep<-transcript_annot[is.na(transcript_annot[[filterName]])==FALSE & transcript_annot[[filterName]]!="",]$ensembl_transcript_id
            to_keep_all<-c(to_keep_all,to_keep)
          }
          else if (filterName=="transcript_appris") {
            #keep transcripts with appropriate transcript_appris values
            if (is.null(appris_options)) {
              #if transcript_appris is given as a filter but appris_options are unspecified, keep principal transcripts only
              message("transcript_appris provided as filter, but no appris_options provided. providing principal appris transcripts (P1-P5)")
              to_keep<-transcript_annot[is.na(transcript_annot$transcript_appris) == FALSE & transcript_annot$transcript_appris !="" & transcript_annot$transcript_appris %in% c("principal1","principal2","principal3","principal4","principal5"),]$ensembl_transcript_id
              to_keep_all<-c(to_keep_all,to_keep)
            }
            else {
              #what to do when appris_options provided.
              #can take both principal (P1-P5) and alternative transcripts (A1-A5).
              #assumes first character of a given option indicates whether it is principal or alternative, so "p" or "P" would indicate principal and "A" or "a" would indicate alternative
              #assumes last character indicates the number (should be between 1 and 5)
              #this way, principal1,Principal1,P1,p1 will all be read as principal1; alternative1,Alternative1, A1,a1 will all be read as alternative1, etc.
              appris_filters_all<-c()
              for (appris_filterName in appris_options) {
                if (substr(appris_filterName,1,1) %in% c("P","p","A","a")==FALSE) {
                  message(paste0("invalid appris filter: ", appris_filterName," - skipping"))
                }
                else {
                  if (substr(appris_filterName,1,1) %in% c("P","p")) {
                    realName<-paste0("principal",substr(appris_filterName,nchar(appris_filterName),nchar(appris_filterName)))
                    if (realName %in% c("principal1","principal2","principal3","principal4","principal5")==FALSE) {
                      message(paste0("invalid appris filter: ", appris_filterName, " - skipping"))
                    }
                    else {
                      appris_filters_all<-c(appris_filters_all,realName)
                    }
                  }
                  else {
                    if (substr(appris_filterName,1,1) %in% c("A","a")) {
                      realName<-paste0("alternative",substr(appris_filterName,nchar(appris_filterName),nchar(appris_filterName)))
                      if (realName %in% c("alternative1", "alternative2","alternative3","alternative4","alternative5")==FALSE) {
                        message(paste0("invalid appris filter: ", appris_filterName, " - skipping"))
                      }
                      else {
                        appris_filters_all<-c(appris_filters_all, realName)
                      }
                    }
                  }
                  appris_filters_all<-unique(appris_filters_all)
                  if (length(appris_filters_all) == 0) {
                    message("no valid appris filters. using principal transcripts")
                    appris_filters_all<-c("principal1", "principal2","principal3","principal4", "principal5")
                  }
                  to_keep<-transcript_annot[is.na(transcript_annot$transcript_appris)==FALSE & transcript_annot$transcript_appris !="" & transcript_annot$transcript_appris %in% appris_filters_all,]$ensembl_transcript_id
                  to_keep_all<-c(to_keep_all,to_keep)
                }
              }
            }
          }
          else if (filterName =="transcript_tsl") {
            if (is.null(tsl_options)) {
              message("no specific transcript support levels provided. will keep all transcripts with any indicated transcript support level")
              to_keep<-transcript_annot[is.na(transcript_annot$transcript_tsl)==FALSE & transcript_annot$transcript_tsl != "",]$ensembl_transcript_id
              to_keep_all<-c(to_keep_all,to_keep)
            }
            else { #use specified transcript support levels (tsl_options) for filtering
              #tsl levels may be given by just numbers, or as tsl1, tsl2, etc
              tsl_all<-c()
              for (tsl in tsl_options) {
                if (nchar(as.character(tsl))>1){ #has characters in addition to support level
                  tsl_real<-tsl
                }
                else { #only the number
                  tsl_real<-paste0("tsl",tsl)
                }
                tsl_all<-c(tsl_all, tsl_real)
              }

              to_keep<-transcript_annot[is.na(transcript_annot$transcript_tsl)==FALSE & transcript_annot$transcript_tsl != "" & substr(transcript_annot$transcript_tsl, 1, 4) %in% tsl_all, ]$ensembl_transcript_id
              to_keep_all<-c(to_keep_all,to_keep)
            }
          }
          transcript_annot<-transcript_annot[transcript_annot$ensembl_transcript_id %in% unique(to_keep_all),]

        }
      }
      else if (canonicalTranscriptOnly==TRUE) {
        transcript_annot<-transcript_annot[is.na(transcript_annot$transcript_is_canonical)==FALSE,]
      }

      #keep only those transcripts that are in range
      transcript_annot$minimum<-pmin(transcript_annot$transcript_start,transcript_annot$transcript_end, na.rm=TRUE)
      transcript_annot<-transcript_annot[transcript_annot$minimum<=max(as.numeric(start),as.numeric(end)),]
      transcript_annot<-transcript_annot[,colnames(transcript_annot) !="minimum"]

      #get exon annotations


      exon_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"),
                                           filters = c("chromosome_name", "start", "end"),
                                           values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                           mart = mart)}, error=function(e) {
                                                   message("getBM failed: ", e$message)
                                                   NULL
                                                 })
      if (is.null(exon_annot)) {
        message("trying useast mirror")
        if(exists("eastMart")==FALSE) {
          eastMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "useast")
        }
        exon_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"),
                                             filters = c("chromosome_name", "start", "end"),
                                             values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                             mart = eastMart)}, error=function(e) {
                                               message("getBM failed: ", e$message)
                                               NULL
                                             })
        if (is.null(exon_annot)) {
          message("trying www mirror")
        }
        if (exists("westMart")==FALSE) {
          westMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "www")
        }

        exon_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"),
                                             filters = c("chromosome_name", "start", "end"),
                                             values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                             mart = westMart)}, error=function(e) {
                                               message("getBM failed: ", e$message)
                                               NULL
                                             })

        if(is.null(exon_annot)) {
          message("trying asia mirror")
          if (exists("asiaMart")==FALSE) {
            asiaMart<-biomaRt::useEnsembl(biomart = "genes",dataset = ensembl_set, mirror = "asia")
          }
          exon_annot<-tryCatch({biomaRt::getBM(attributes = c("ensembl_transcript_id", "exon_chrom_start", "exon_chrom_end"),
                                               filters = c("chromosome_name", "start", "end"),
                                               values = list(chr, min(as.numeric(start),as.numeric(end)),max(as.numeric(start),as.numeric(end))),
                                               mart = asiaMart)}, error=function(e) {
                                                 message("getBM failed: ", e$message)
                                                 NULL
                                               })

          if(is.null(exon_annot)) {
            stop("failed to get exon annotations")
          }

        }
      }





      exon_annot<-exon_annot[exon_annot$ensembl_transcript_id %in% transcript_annot$ensembl_transcript_id,]


      names(exon_annot)[names(exon_annot)=="ensembl_transcript_id"]<-"transcript_id"
      names(exon_annot)[names(exon_annot)=="exon_chrom_start"]<-"start"
      names(exon_annot)[names(exon_annot)=="exon_chrom_end"]<-"end"

      names(transcript_annot)[names(transcript_annot)=="ensembl_transcript_id"]<-"transcript_id"
      names(transcript_annot)[names(transcript_annot)=="transcript_start"]<-"start"
      names(transcript_annot)[names(transcript_annot)=="transcript_end"]<-"end"



    }

    #getting relevant transcripts from custom annotation
    else {
      #get transcripts in region of interest
      transcript_gr<-annotation[annotation$type=="transcript"]
      overlap_transcript_gr <- GenomicRanges::findOverlaps(transcript_gr, region_gr)
      transcriptsInRegion <- transcript_gr[S4Vectors::queryHits(overlap_transcript_gr)]
      transcript_annot<-as.data.frame(transcriptsInRegion)

      if ("transcript_id" %in% colnames(transcript_annot) == FALSE) {
        stop("error: 'transcript_id' not found in provided annotation")
      }

      gene_annot_state<-FALSE
      if (base::nrow(gene_annot) !=0) {
        gene_annot_state<-TRUE
        transcript_annot<-transcript_annot[transcript_annot[[gene_symbol]] %in% gene_annot[[gene_symbol]] | transcript_annot$transcript_id %in% transcript_list ,]
      }
      else {
        transcript_annot<-transcript_annot[transcript_annot$transcript_id %in% transcript_list,]
      }

      #keep relevant transcripts
      if (is.null(transcript_list)==FALSE) {
        transcript_annot<-transcript_annot[transcript_annot$transcript_id %in% transcript_list,]
        if (base::nrow(transcript_annot)==0) {
          message("no transcripts found from provided list. will not display any transcript info")
          transcript_annot<-NULL
        }
        else if (base::nrow(transcript_annot)<length(transcript_list)) {
          message(paste(c("transcripts not found:",transcript_list[transcript_list %in% transcript_annot$transcript_id == FALSE], "- ignoring" ), collapse=" "))
        }
      }
      else if (supportedTranscriptsOnly==TRUE & canonicalTranscriptOnly==FALSE) {
        #only use available filters
        transcript_filters_present<-transcript_filters[transcript_filters %in% colnames(transcript_annot)]
        if (length(transcript_filters_present) == 0) {
          if (is.null(tsl_options) & is.null(tag_options)) {
            message("no transcript filters found. using all transcripts in annotation.")
          }
        }

        to_keep_all<-c()
        for (filterName in transcript_filters_present) {
          if (filterName != "tag" & filterName != "transcript_support_level") {
            #for these filters, take anything non-empty
            to_keep<-transcript_annot[is.na(transcript_annot[[filterName]]) == FALSE & transcript_annot[[filterName]] %in% c("NA","")==FALSE,]$transcript_id
            to_keep_all<-c(to_keep_all, to_keep)
          }
          else if (filterName == "tag") {
            if (is.null(tag_options)) {
              message("filter 'tag' found but no tag_options are specified. only removing transcripts without any tag")
              to_keep<-transcript_annot[is.na(transcript_annot$tag)==FALSE & transcript_annot$tag %in% c("NA","")==FALSE ,]$transcript_id
            }
            else { #when tag_options specified
              to_keep<-transcript_annot[is.na(transcript_annot$tag)==FALSE & transcript_annot$tag %in% c("NA","")==FALSE & transcript_annot$tag %in% tag_options,]$transcript_id

            }
            to_keep_all<-c(to_keep_all, to_keep)
          }
          else if (filterName=="transcript_support_level") {
            if (is.null(tsl_options)) {
              message("no specific transcript support levels provided. will keep all transcripts with any indicated transcript support level")
              to_keep<-transcript_annot[is.na(transcript_annot$transcript_support_level)==FALSE & transcript_annot$transcript_support_level %in% c("NA","") ==FALSE,]$ensembl_transcript_id
              to_keep_all<-c(to_keep_all,to_keep)
            }
            else { #when tsl_options specified
              to_keep<-transcript_annot[is.na(transcript_annot$transcript_support_level)==FALSE & transcript_annot$transcript_support_level %in% c("NA","")==FALSE &
                                          (transcript_annot$transcript_support_level %in% tsl_options | transcript_annot$transcript_support_level %in% paste0("tsl",tsl_options) | paste0(tsl,transcript_annot$transcript_support_level) %in% tsl_options),]$transcript_id
              to_keep_all<-c(to_keep_all,to_keep)

            }

          }

        }

      }

      else if (canonicalTranscriptOnly==TRUE) {
        #will first search for column names in transcript_annot that would indicate information about canonical transcript. if it finds none, it will see if there are any transcripts where tag=MANE_select.
        to_keep_all<-c()
        potential_canonical_names<-c("transcript_is_canonical", "canonical_transcript", "is_canonical","canonical", "MANE_select", "Ensembl_canonical","EnsEMBL_canonical")
        allValues<-colnames(transcript_annot)
        if (length(allValues[allValues %in% potential_canonical_names])==0 ) { #what to do if none of the canonical names are present - use MANE_select transcript if available as a tag
          if ("tag" %in% allValues) {
            if (base::nrow(transcript_annot[transcript_annot$tag=="MANE_Select",])>0) {
              message("using MANE_select as canonical transcript")
              to_keep<-transcript_annot[is.na(transcript_annot$tag)==FALSE & transcript_annot$tag=="MANE_select",]$transcript_id
              to_keep_all<-c(to_keep_all, to_keep)
            }
            else { #no mane_select transcripts
              message("no information about canonical transcript found. will not display any transcripts")
            }
          }
          else { #tag filter unavailable
            message("no information about canonical transcript found. will not display any transcripts")
          }
        }
        else {
          canonicalFilters<-potential_canonical_names[potential_canonical_names %in% allValues]
          for (canFilt in canonicalFilters) {
            to_keep<-transcript_annot[is.na(transcript_annot[[canFilt]])==FALSE & transcript_annot[[canFilt]] %in% c(FALSE, "NA","",0,"0")==FALSE,]$transcript_id
            to_keep_all<-c(to_keep_all, to_keep)
          }
        }
        transcript_annot<-transcript_annot[transcript_annot$transcript_id %in% to_keep_all,]
      }

      #get exons in region of interest
      exon_gr<-annotation[annotation$type=="exon"]
      overlap_exon_gr <- GenomicRanges::findOverlaps(exon_gr, region_gr)
      exonsInRegion <- exon_gr[S4Vectors::queryHits(overlap_exon_gr)]
      exon_annot<-as.data.frame(exonsInRegion)

      names(exon_annot)[names(exon_annot)=="seqnames"]<-"chromosome_name"
      names(transcript_annot)[names(transcript_annot)=="seqnames"]<-"chromosome_name"

    }

    #keep coordinates within range
    transcript_annot$minimum<-pmin(transcript_annot$start,transcript_annot$end, na.rm=TRUE)

    transcript_annot<-transcript_annot[transcript_annot$minimum<=max(as.numeric(start),as.numeric(end)),]

    transcript_annot<-transcript_annot[,colnames(transcript_annot) !="minimum"]


    if (base::nrow(transcript_annot) !=0) {
      transcript_annot$rowNum<-1:base::nrow(transcript_annot)
      for(transcriptNum in transcript_annot$rowNum) {
        #plotting transcripts with start/end outside of range may result in error. adjust transcript coordinates to fall in range, if necessary:
        startCoord<-transcript_annot[transcript_annot$rowNum==geneNum,]$start
        endCoord<-transcript_annot[transcript_annot$rowNum==geneNum,]$end
        transcript_annot[transcript_annot$rowNum==geneNum,]$start<-max(startCoord,min(as.numeric(start),as.numeric(end)))
        transcript_annot[transcript_annot$rowNum==geneNum,]$end<-min(endCoord,max(as.numeric(start),as.numeric(end)))
      }

      if (base::nrow(exon_annot) !=0) {
        exon_annot$rowNum<-1:base::nrow(exon_annot)

      }

    }

    if (gene_annot_state==TRUE) {
      #put it all together for plotting
      allGenes<-unique(gene_annot[[gene_symbol]])

      graphLoc<-0

      allInfo<-data.frame()
      for (gene in allGenes) {
        geneStart<-gene_annot[gene_annot[[gene_symbol]]==gene,]$start
        geneEnd<-gene_annot[gene_annot[[gene_symbol]]==gene,]$end
        strand<-gene_annot[gene_annot[[gene_symbol]]==gene,]$strand
        geneScore<-graphLoc

        gene_df<-data.frame(chromosome_name=chr, start=geneStart, end=geneEnd, strand=strand,
                            type="Gene",graph_location=geneScore, geneName=gene,transcriptName=NA, label=gene,
                            start_real=gene_annot[gene_annot[[gene_symbol]]==gene,]$start_real,
                            end_real=gene_annot[gene_annot[[gene_symbol]]==gene,]$end_real)
        allInfo<-rbind(allInfo, gene_df)


        #get transcripts
        transcript_gene<-transcript_annot[transcript_annot[[gene_symbol]]==gene,]

        if (base::nrow(transcript_gene)>0) {
          gene_annot$graph_location<-seq(0, by=-3, length.out=base::nrow(gene_annot))

          transcript_gene$graph_location<-seq(geneScore-2, by=-2, length.out=base::nrow(transcript_gene))
          graphLoc<-min(transcript_gene$graph_location)-3

          transcript_df<-data.frame(chromosome_name=chr, start=transcript_gene$start, end=transcript_gene$end,
                                    strand=strand, type="Transcript", graph_location=transcript_gene$graph_location, geneName=gene,
                                    transcriptName=transcript_gene$transcript_id, label=transcript_gene$transcript_id,
                                    start_real=transcript_gene$start, end_real=transcript_gene$end)
          transcript_df$start[transcript_df$start_real<min(start,end)]<-min(start,end)
          transcript_df$end[transcript_df$end_real>max(start,end)]<-max(start,end)
          allInfo<-rbind(allInfo,transcript_df)
        }

        #get exons
        for (transcript_name in transcript_gene$transcript_id ) {
          exon_transcript<-exon_annot[exon_annot$transcript_id == transcript_name,]
          exon_transcript$graph_location<-transcript_gene[transcript_gene$transcript_id==transcript_name,]$graph_location

          exon_df<-data.frame(chromosome_name=chr, start=exon_transcript$start, end=exon_transcript$end,
                              strand=strand, type="Exon", graph_location=exon_transcript$graph_location, geneName=gene,
                              transcriptName=transcript_name,label=NA,
                              start_real=exon_transcript$start, end_real=exon_transcript$end)
          exon_df$start[exon_df$start_real<min(start,end)]<-min(start,end)
          exon_df$end[exon_df$end_real>max(start,end)]<-max(start,end)
          allInfo<-rbind(allInfo, exon_df)
        }

      }


      allInfo<-assign_y_levels_transcripts(allInfo)

      allInfo$start<-allInfo$start/1000000
      allInfo$end<-allInfo$end/1000000

      allInfo$midX<-(allInfo$start+allInfo$end)/2
      allInfo$txt_label_loc<-allInfo$graph_location-0.75
      allInfo$trackType<-"gene"
      allInfo$plot_start<-min(start,end)
      allInfo$plot_end<-max(start,end)
      allInfo$includeTxtNames<-includeTxtNames
      arrow_properties_gene <- grid::arrow(type = "closed", angle = 20, length = grid::unit(0.3, "cm"))


      p1<-ggplot2::ggplot(allInfo, ggplot2::aes(x = .data$start, y = .data$graph_location)) +
        ggplot2::geom_segment(data=allInfo[allInfo$type=="Transcript",], ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                     linewidth=0.5, lineend = "butt", linejoin = "mitre")+
        ggplot2::geom_segment(data=allInfo[allInfo$type=="Exon",], ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                     linewidth=1.5, lineend = "butt", linejoin = "mitre")+
        ggplot2::geom_segment(data=allInfo[allInfo$type=="Gene",],ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                     linewidth=3, arrow=arrow_properties_gene, lineend = "butt", linejoin = "mitre")+
        ggplot2::geom_label(data=allInfo[allInfo$type=="Gene",],color="black", fill="white",ggplot2::aes(x=.data$midX, y=.data$graph_location, label=.data$label),
                   fontface = "italic",size = fontSize_text)+
        ggplot2::xlim(min(start,end)/1000000, max(start,end)/1000000)+
        ggplot2::ylim(min(allInfo$graph_location)-1,2)+ggplot2::theme_classic()+ggplot2::ggtitle(graphTitle)+
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              axis.line.y = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5))+
        ggplot2::xlab("Location (Mb)")

      if (includeTxtNames==TRUE) {
        p1<-p1+ggplot2::geom_text(data=allInfo[allInfo$type=="Transcript",], color="black", ggplot2::aes(x=.data$midX,y=.data$txt_label_loc, label=.data$label), size=fontSize_text)
      }


    }

    else {
      gene_annot<-data.frame(chromosome_name=chr, start=as.numeric(start)/1000000, end=as.numeric(end)/1000000, gene_id=NA, strand=NA,
                             type=NA, label=NA, trackColor=NA, rowNum=NA, start_real=start, end_real=end, graph_location=0, midX=((as.numeric(start)+as.numeric(end))/2)/1000000, trackType="gene", plot_start=start, plot_end=end )
      p1<-ggplot2::ggplot(gene_annot)+ggplot2::geom_blank() + ggplot2::xlim(min(as.numeric(start),as.numeric(end))/1000000, max(as.numeric(start),as.numeric(end))/1000000) +
        ggplot2::ylim(1,2) +ggplot2::theme_classic()+ ggplot2::ggtitle(graphTitle)+
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              axis.line.y = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5))+
        ggplot2::xlab("Location (Mb)")

    }
  }

  return(list(figure=p1, genePlot=p1))

}



#' plot coverage tracks
#' @param genomicLoc character vector of the form c(chromosome, start coordinate, end coordinate) (ex: c("chrX", 24040226, 24232664))
#' @param covFiles character vector containing a list of file paths to any coverage files (typically in bigWig or bedGraph format)
#' @param covTrackNames character vector containing the label of each coverage track. Should be the same length as covFiles and in the same order. By default it will label them as "Cov_1", "Cov_2", etc
#' @param covTrackColors character vector containing the color of each coverage track. If only one color is provided, all tracks will be that color. If different tracks must be different colors, specify the color for each track in order. Default is "black"
#' @param ymin numeric, minimum value for y-axis for coverage plot
#' @param ymax numeric, maximum value for y-axis for coverage plot
#' @param fillArea boolean, set this to FALSE if you do not want to fill in the area under your coverage tracks. Default is TRUE
#' @param logScale boolean, set this to TRUE if you want data plotted on a log scale. Default is FALSE
#' @param rasterizePlot boolean, set this to TRUE if you want to rasterize the coverage plot. This can be useful as coverage files can be rather large, and that can mess up the plot when you save it as a vector file. Default is FALSE
#' @param fontSize a numeric for desired font size. Default is 9.
#' @return ggplot of coverage tracks

plot_coverage<-function(genomicLoc, covFiles, covTrackNames=NULL,covTrackColors="black",ymin=NULL, ymax=NULL, fillArea=TRUE,logScale=FALSE, rasterizePlot=FALSE, fontSize=9) {
  if (length(genomicLoc) != 3) {
    stop("error in plot_coverage: please provide genomic coordinates")
  }
  #set up covFiles_df
  if (is.null(covTrackNames)) {
    covTrackNames<-paste0("Cov_", 1:length(covFiles)) #will name the tracks in the order they came
  }
  else {
    if (length(covTrackNames) != length(covFiles)) {
      stop("error: please make sure the number of track names provided is the same as the number of tracks")
    }
  }

  if (length(covTrackColors) != 1 & length(covTrackColors) !=length(covFiles)) {
    stop("error: please either select one color for all tracks, or provide the same number of colors as number of tracks")
  }
  covFiles_df<-data.frame(fileNames=covFiles, displayNames=covTrackNames, colorNames=covTrackColors)
  start<-as.numeric(genomicLoc[2])
  end<-as.numeric(genomicLoc[3])
  #the chromosome name might include "chr" or not -- make sure that this is handled
  if (grepl("chr", genomicLoc[1])==TRUE) {
    graphTitle<-paste0(genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
    gene_info<-data.frame(chromosome_name=genomicLoc[1], start_position=genomicLoc[2], end_position=genomicLoc[3])
    chr <- substr(genomicLoc[1], 4, nchar(genomicLoc[1]))

  }
  else {
    if (suppressWarnings(is.na(as.numeric(genomicLoc[1]))==FALSE) | genomicLoc[1] %in% c("X","Y","Z","W")) {
      graphTitle<-paste0("chr",genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
      gene_info<-data.frame(chromosome_name=paste0("chr",genomicLoc[1]), start_position=genomicLoc[2], end_position=genomicLoc[3])
      chr<-genomicLoc[1]
    }

    else {
      graphTitle<-paste0(genomicLoc[1],":",genomicLoc[2],"-",genomicLoc[3]) #used for the title of the plot
      gene_info<-data.frame(chromosome_name=genomicLoc[1], start_position=genomicLoc[2], end_position=genomicLoc[3])
      chr<-as.character(genomicLoc[1])
    }

  }


  # Import trackfiles
  allSamples<-data.frame()
  for (trackFile in covFiles_df$fileNames) {
    chromosomes<-get_chrom_names(trackFile)
    if (chr %in% chromosomes) {
      gene_range<-GenomicRanges::GRanges(seqnames = chr,ranges = IRanges::IRanges(start = as.numeric(start),end = as.numeric(end)))
    }
    else {
      if (paste0("chr",chr) %in% chromosomes) {
        gene_range<-GenomicRanges::GRanges(seqnames=paste0("chr",chr), ranges=IRanges::IRanges(start=as.numeric(start), end=as.numeric(end)))
      }
      else {
        stop("chromosome not found in coverage file")
      }
    }

    if (tools::file_ext(trackFile) %in% c("bdg","bg", "bedgraph", "bedGraph")) {
      bigwig_df<-read_bedgraph_region(trackFile, as.character(GenomeInfoDb::seqnames(gene_range)), as.numeric(start), as.numeric(end))

    }
    else {
      bigwig_gRanges<-rtracklayer::import(trackFile, which=gene_range)
      bigwig_df<-as.data.frame(bigwig_gRanges)
    }

    bigwig_df$colorNames<-covFiles_df[covFiles_df$fileNames==trackFile,]$colorNames
    bigwig_df$displayNames<-covFiles_df[covFiles_df$fileNames==trackFile,]$displayNames
    allSamples<-rbind(allSamples, bigwig_df)
  }

  allSamples<-data.frame(allSamples, stringsAsFactors = FALSE, check.names = FALSE)

  allSamples$start<-as.numeric(allSamples$start)/1000000
  allSamples$end<-as.numeric(allSamples$end)/1000000
  names(allSamples)[names(allSamples)=="seqnames"]<-"chr"

  allSamples$logTransform<-logScale
  allSamples$rasterize<-rasterizePlot
  allSamples$fillArea<-fillArea


  allSamples$displayNames<-factor(allSamples$displayNames, levels=covTrackNames, ordered=TRUE)

  if (is.null(ymin)) {
    ymin<-min(allSamples$score)
  }
  if (is.null(ymax)) {
    ymax<-max(allSamples$score)
  }

  #if log transforming
  if (logScale==TRUE) {
    #need to drop all data points less than or equal to 0
    allSamples<-allSamples[allSamples$score >0 ,]
    #plotting range also can't have values less than or equal to zero

    if (ymin <= 0) {
      if (ymax > 1) { #set minimum value to 1 as long as ymax is bigger than 1
        message("invalid ymin value for log scale. setting to 1")
        ymin<-1
      }
      else {
        if (ymax > 0) { #if ymax is between 0 and 1, set ymin to be ymax/10
          ymin<-ymax/10
          message(paste(c("invalid ymin value for log scale. setting to", ymin), sep=" "))
        }
        else { #if ymax is also negative
          stop("error: invalid y axis range for log scale")
        }
      }

    }
    #plot the log transform
    if (fillArea==TRUE) {
      p1<-ggplot2::ggplot(allSamples, ggplot2::aes(x = .data$start, y = .data$score, color=.data$colorNames, fill=.data$colorNames)) +
        ggplot2::geom_line() + ggplot2::geom_area()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
        ggplot2::facet_grid(displayNames ~ ., switch = "y", scales = "free_y", drop=FALSE)+ggplot2::ggtitle(graphTitle)+ggplot2::xlab("Location (Mb)")+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),strip.text = ggplot2::element_text(size = fontSize), axis.text.y=ggplot2::element_text(size=fontSize))+
        ggplot2::scale_y_continuous(trans='log2',n.breaks = 3)+
        ggplot2::coord_cartesian(xlim = c(min(as.numeric(start),as.numeric(end))/1000000, max(as.numeric(start),as.numeric(end))/1000000), ylim = c(ymin, ymax))+ggplot2::theme(axis.title.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=fontSize))
    }
    else {
      p1<-ggplot2::ggplot(allSamples, ggplot2::aes(x = .data$start, xend=.data$end, y = .data$score, color=.data$colorNames, fill=.data$colorNames)) +
        ggplot2::geom_segment()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
        ggplot2::facet_grid(displayNames ~ ., switch = "y", scales = "free_y", drop=FALSE)+ggplot2::ggtitle(graphTitle)+ggplot2::xlab("Location (Mb)")+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),strip.text = ggplot2::element_text(size = fontSize), axis.text.y=ggplot2::element_text(size=fontSize))+
        ggplot2::scale_y_continuous(trans='log2',n.breaks = 3)+
        ggplot2::coord_cartesian(xlim = c(min(as.numeric(start),as.numeric(end))/1000000, max(as.numeric(start),as.numeric(end))/1000000), ylim = c(ymin, ymax))+ggplot2::theme(axis.title.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=fontSize))
    }

  }

  else { #plot without log transform
    if (fillArea==TRUE) {


      start_num <- min(as.numeric(unclass(start)), as.numeric(unclass(end)))
      end_num   <- max(as.numeric(unclass(start)), as.numeric(unclass(end)))

      p1<-ggplot2::ggplot(allSamples, ggplot2::aes(x = .data$start, y = .data$score, color=.data$colorNames, fill=.data$colorNames,group = .data$colorNames)) +
        ggplot2::geom_line() + ggplot2::geom_area()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
        ggplot2::facet_grid(displayNames ~ ., switch = "y", scales = "free_y")+
        ggplot2::ggtitle(graphTitle)+ggplot2::xlab("Location (Mb)")+
        ggplot2::scale_y_continuous(n.breaks = 3)+ggplot2::coord_cartesian(xlim = c(start_num/1000000, end_num/1000000), ylim = c(ymin, ymax))+
        ggplot2::theme(axis.title.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=fontSize),axis.text.x=ggplot2::element_text(size=fontSize),plot.title = ggplot2::element_text(hjust = 0.5),strip.text = ggplot2::element_text(size = fontSize))
    }
    else {
      p1<-ggplot2::ggplot(allSamples, ggplot2::aes(x = .data$start, xend=.data$end, y = .data$score, color=.data$colorNames, fill=.data$colorNames)) +
        ggplot2::geom_segment()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
        ggplot2::facet_grid(displayNames ~ ., switch = "y", scales = "free_y", drop=FALSE)+ggplot2::ggtitle(graphTitle)+ggplot2::xlab("Location (Mb)")+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),strip.text = ggplot2::element_text(size = fontSize),axis.text.y=ggplot2::element_text(size=fontSize))+ggplot2::scale_y_continuous(n.breaks = 3)+
        ggplot2::coord_cartesian(xlim = c(min(as.numeric(start),as.numeric(end))/1000000, max(as.numeric(start),as.numeric(end))/1000000), ylim = c(ymin, ymax))+ggplot2::theme(axis.title.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=fontSize))
    }
  }
  if (rasterizePlot==TRUE) {
    rast_layers <- lapply(p1$layers, function(layer) {
      if (inherits(layer$geom, "GeomArea") || inherits(layer$geom, "GeomLine") || inherits(layer$geom,"GeomSegment")) {
        ggrastr::rasterise(layer, dpi = 300)
      } else {
        layer
      }
    })

    p1$layers<-rast_layers

  }


  return(p1)
}


#' plot peak tracks
#' @param genomicLoc character vector of the form c(chromosome, start coordinate, end coordinate) (ex: c("chrX", 24040226, 24232664))
#' @param peakFiles character vector containing a list of file paths to any peak files (.bed format)
#' @param peakTrackNames character vector containing the label of each peak track. Should be the same length as 'peakFiles' and in the same order. By default it will label them as "Peaks_1", "Peaks_2", etc
#' @param peakTrackColors character vector containing color of each peak track. If only one color is provided, all tracks will be that color. If different tracks must be different colors, specify the color for each track in order. Default is "black"
#' @param labelAllPeaks boolean, set to TRUE if each individual peak should be named. assumes names are provided in the fourth column of a given .bed file
#' @param specialPeaks character vector listing the names of specific peaks you want to either label or put in a different color
#' @param labelSpecialPeaks boolean, set to TRUE if each special peak should be named. requires 'specialPeaks' to be specified
#' @param specialPeakColors character vector containing the color(s) you want to give your special peaks. requires 'specialPeaks' to be specified
#' @param labelStrand boolean, set to TRUE if the strand of each peak should be indicated. assumes strand information is in the sixth column of a given .bed file
#' @param strandColors character vector of length 2. First is the color for the + strand, the second is the color for the - strand. overrides 'peakTrackColors' but does NOT override 'specialPeakColors'
#' @param fontSize a numeric for desired font size. Default is 9.
#' @return ggplot of peak tracks

plot_peaks<-function(genomicLoc, peakFiles, peakTrackNames=NULL,
                     labelAllPeaks=FALSE, peakTrackColors="black",
                     specialPeaks=NULL, labelSpecialPeaks=FALSE, specialPeakColors=NULL,
                     labelStrand=FALSE, strandColors=NULL, fontSize=9) {
  fontSize_text<-fontSize/ggplot2::.pt
  if (length(genomicLoc) != 3) {
    stop("error in plot_peaks: please provide genomic coordinates")
  }

  #if names for tracks not provided, use "Track_1", "Track_2", etc.
  if (is.null(peakTrackNames)) {
    numberOfTracks<-1:length(peakFiles)
    peakTrackNames<-paste0("Peaks_", numberOfTracks)
  }

  else {
    #number of names for peak tracks should be the same as the number of peak files provided
    if (length(peakTrackNames) != length(peakFiles)) {
      stop("error: number of peak files is not equal to the number of track names")
    }
  }

  if (length(peakTrackColors) !=1 & length(peakTrackColors) !=length(peakFiles)) { #either one color for all tracks, or a color specified for every individual track
    stop("error: number of peak files is not equal to the number of peak colors" )
  }

  chr<-genomicLoc[1]
  if (grepl("chr",chr)==TRUE) {
    chr<-substr(chr, 4, nchar(chr))
  }
  start<-as.numeric(genomicLoc[2])
  end<-as.numeric(genomicLoc[3])

  graphTitle<-paste0("chr",chr,":",start,"-",end)


  allPeaks<-data.frame()
  for (samplePeaks in peakFiles) {
    if (is.data.frame(samplePeaks)==FALSE) {
      #check if file exists
      if (file.exists(samplePeaks) == FALSE){
        stop(paste0("error: ", samplePeaks, " peak file not found"))

      }
      else {
        peaks<-read.delim(samplePeaks, check.names=FALSE, stringsAsFactors = FALSE, header = FALSE)

      }
    }
    else { #when the provided object is a dataframe rather than a file name
      peaks<-samplePeaks
    }

    names(peaks)[1:3] <- c("chr", "start", "end")
    if (ncol(peaks)>3) { #if peak names or strands are provided, include that information
      names(peaks)[4]<-"name"
      if (ncol(peaks)>5) {
        names(peaks)[6]<-"strand"
      }
    }
    peaks<-peaks[,colnames(peaks) %in% c("chr", "start","end","name","strand")]
    if ("name" %in% colnames(peaks)==FALSE) {
      peaks$name<-""
    }
    if ("strand" %in% colnames(peaks)==FALSE) {
      peaks$strand<-""
    }
    #keep only peaks in the region
    peaks<-peaks[peaks$chr==chr | peaks$chr==paste0("chr",chr),]
    peaks<-peaks[peaks$end>as.numeric(start) & peaks$start<as.numeric(end),]


    if (nrow(peaks)!=0) {
      peaks$displayNames<-NA
      #add track name info
      fileIndex<-which(sapply(peakFiles, `==`, samplePeaks))
      peaks$displayNames<-peakTrackNames[fileIndex]
      #add Color info
      peaks$colorNames<-NA
      if (length(peakTrackColors) == 1) {
        peaks$colorNames<-peakTrackColors
      }
      else {
        peaks$colorNames<-peakTrackColors[fileIndex]
      }
      allPeaks<-rbind(allPeaks,peaks)

    }

  }


  if (nrow(allPeaks)==0) {
    message("no peaks in region")
    p1<-NULL
  }

  else {
    #coloring by strand
    if (is.null(strandColors)==FALSE) {
      if (length(strandColors)!=2) {
        message("error: strandColors should have two different colors - the first for the + strand, the second for the - strand. this does not seem to be the case. peaks will not be colored according to strand")
      }
      else {
        allPeaks$colorNames[is.na(allPeaks$strand)==FALSE &(allPeaks$strand=="+" | allPeaks$strand == 1)]<-strandColors[1]
        allPeaks$colorNames[is.na(allPeaks$strand)==FALSE & (allPeaks$strand=="-" | allPeaks$strand==(-1))]<-strandColors[2]

      }
    }

    #coloring special peaks
    if (is.null(specialPeakColors)==FALSE) {
      if (length(specialPeakColors) !=1 & length(specialPeakColors) !=length(specialPeaks)) {
        message("error: number of specialPeakColors does not match number of specialPeaks. will not color any specialPeaks")
      }
      else {
        if (length(specialPeakColors)==1) {
          allPeaks$colorNames[allPeaks$name %in% specialPeaks]<-specialPeakColors
        }
        else {
          for (specialPeakName in specialPeaks) {
            peakColorIndex<-which(sapply(specialPeaks, `==`, specialPeakName))
            peakColorName<-specialPeakColors[peakColorIndex]
            allPeaks$colorNames[allPeaks$name==specialPeakName]<-peakColorName
          }
        }
      }
    }

    #add labels, if necessary
    allPeaks$label<-NA
    if (labelAllPeaks == TRUE) {
      allPeaks$label<-allPeaks$name
    }

    if (labelSpecialPeaks == TRUE) {
      allPeaks$label[allPeaks$name %in% specialPeaks]<-allPeaks[allPeaks$name %in% specialPeaks,]$name
    }

    allPeaks$strandLabel<-NA
    if (labelStrand==TRUE) {
      allPeaks$strandLabel<-allPeaks$strand
    }

    allPeaks$trackType<-"Peaks"

    allPeaks$start<-allPeaks$start/1000000
    allPeaks$end<-allPeaks$end/1000000

    #labels will be put in the middle of peaks. however, if part of the peak is out of range, that might get messed up. adjust accordingly
    allPeaks$start_forLabel<-allPeaks$start
    allPeaks$start_forLabel[allPeaks$start_forLabel<as.numeric(start)/1000000]<-as.numeric(start)/1000000

    allPeaks$end_forLabel<-allPeaks$end
    allPeaks$end_forLabel[allPeaks$end_forLabel>as.numeric(end)/1000000]<-as.numeric(end)/1000000

    allPeaks$midX<-(allPeaks$start_forLabel + allPeaks$end_forLabel)/2



    allPeaks$graph_location<-1


    #build the plot
    allPeaks$displayNames<-as.character(allPeaks$displayNames)
    allPeaks$displayNames<-factor(allPeaks$displayNames, levels=peakTrackNames)
    allPeaks$labelStrand<-labelStrand
    allPeaks$labelAllPeaks<-labelAllPeaks
    allPeaks$labelSpecialPeaks<-labelSpecialPeaks

    p1<-ggplot2::ggplot(allPeaks, ggplot2::aes(x = .data$start, color=.data$colorNames)) +
      ggplot2::geom_segment(data=allPeaks, ggplot2::aes(x=.data$start, xend=.data$end, y=1, yend=1,
                                      color=.data$colorNames),linewidth=3, lineend = "butt", linejoin = "mitre")+
      ggplot2::scale_color_identity()+ggplot2::theme_classic()+
      ggplot2::coord_cartesian(xlim = c(min(as.numeric(start),as.numeric(end))/1000000, max(as.numeric(start),as.numeric(end))/1000000), ylim=c(0,2))+
      ggplot2::facet_grid(displayNames ~ ., switch = "y", scales = "free_y", drop=FALSE)+
      ggplot2::ggtitle(graphTitle)+ggplot2::xlab("Location (Mb)")+
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5),strip.text = ggplot2::element_text(size = fontSize))

    if (labelStrand == TRUE) {

      p1<-p1+ggplot2::geom_text(data=allPeaks, ggplot2::aes(x=.data$midX, label=.data$strand, y=.data$graph_location), size=fontSize_text, color="black", vjust = 0.5)

    }

    if (labelAllPeaks == TRUE | labelSpecialPeaks== TRUE) {

      p1<-p1+ggrepel::geom_text_repel(data=allPeaks, y=0.85, ggplot2::aes(x=.data$midX, label=.data$label), size=fontSize_text, color="black", nudge_y=-0.2, direction="y", segment.size=0.5,na.rm = TRUE)


    }

  }



  return(p1)

}


#' plot loop tracks
#' @param genomicLoc character vector of the form c(chromosome, start coordinate, end coordinate) (ex: c("chrX", 24040226, 24232664))
#' @param loopFiles character vector containing a list of file paths to any contact/loop files (.bedpe format)
#' @param loopTrackNames character vector containing the label of each loop track. Should be the same length as 'contactFiles' and in the same order. By default it will label them as "Loops_1", "Loops_2", etc
#' @param loopTrackColors character vector containing color of each loop track. If only one color is provided, all tracks will be that color. If different tracks must be different colors, specify the color for each track in order. Default is 'black'
#' @param lineSize numeric, width of lines used to draw loops. Default is 0.8
#' @param alpha numeric, alpha value (transparency) for loops. Default is 0.8
#' @param minScore numeric; any loops with scores lower than this will be thinner and more transparent, making the loops that pass this threshold more visible
#' @param rasterizePlot boolean, set this to TRUE if you want to rasterize the loop plot. This can be useful as loop files can be rather large, and that can mess up the plot when you save it as a vector file. Default is FALSE
#' @param specialLoops character vector listing the names of specific loops you want to put in a different color
#' @param specialLoopColors character vector containing the color(s) you want to give your special loops. requires 'specialLoops' to be specified
#' @param loop_orientation either "above" or "below". "above" will draw loops above a horizontal axis; "below" will draw loops below a horizontal axis. Default is "above'
#' @param fontSize a numeric for desired font size. Default is 9.
#' @return ggplot of loop tracks

plot_loops<-function(genomicLoc, loopFiles, loopTrackNames=NULL,
                     loopTrackColors="black", lineSize=0.8, alpha=0.8,rasterizePlot=FALSE,
                     specialLoops=NULL, specialLoopColors=NULL, minScore=NULL, loop_orientation="above", fontSize=9) {
  if (length(genomicLoc) != 3) {
    stop("error in plot_loops: please provide genomic coordinates")
  }

  if (is.null(loopTrackNames)==FALSE & length(loopTrackNames) != length(loopFiles)) {
    stop("error: number of loop files is not equal to the number of loop names")
  }

  if (is.null(loopTrackNames)) {
    numberOfTracks<-1:length(loopFiles)
    loopTrackNames<-paste0("Loops_", numberOfTracks)
  }

  else {
    if (length(loopTrackNames) != length(loopFiles)) {
      stop("error: number of loop files is not equal to the number of names provided")
    }
  }

  if (length(loopTrackColors) != 1 & length(loopTrackColors) !=length(loopFiles)) {
    stop("error: number of loop files is not equal to the number of loop colors")
  }

  chr<-genomicLoc[1]
  if (grepl("chr",chr)==TRUE) {
    chr<-substr(chr, 4, nchar(chr))
  }
  start<-as.numeric(genomicLoc[2])
  end<-as.numeric(genomicLoc[3])

  graphTitle<-paste0("chr",chr,":",start,"-",end)

  #read in loop files
  allLoops<-data.frame()
  for (sampleLoops in loopFiles) {

    if (is.data.frame(sampleLoops)==FALSE) {
      if (file.exists(sampleLoops)==FALSE) {
        message(paste0("error: ", sampleLoops, " loop file not found. skipping ..."))
      }
      else {
        loops<-read.delim(sampleLoops, check.names=FALSE, stringsAsFactors = FALSE, header=FALSE)
        names(loops)[1:6]<-c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
      }
    }
    else {
      loops<-sampleLoops
    }

    if (ncol(loops)>6) {
      if (is.numeric(loops[,7])) {
        names(loops)[7]<-"score"
      }
      else {
        names(loops)[7]<-"name"
      }

      if (ncol(loops)>7) {
        if (is.numeric(loops[,8])) {
          names(loops)[8]<-"score"
        }
        else {
          names(loops)[8]<-"name"
        }
      }
    }

    if ("score" %in% colnames(loops)==FALSE) {
      loops$score<-NA
    }
    if ("name" %in% colnames(loops)==FALSE) {
      loops$name<-""
    }
    loops<-loops[,c("chrom1","start1","end1","chrom2","start2","end2","name","score")]

    #keep only loops on the chromosome of interest
    loops<-loops[(loops$chrom1 == chr | loops$chrom1 ==paste0("chr",chr)) &
                   (loops$chrom2 == chr | loops$chrom2 ==paste0("chr",chr)),]
    loops$displayNames<-NA

    #add track name info
    fileIndex<-which(sapply(loopFiles, `==`, sampleLoops))
    loops$displayNames<-loopTrackNames[fileIndex]

    #add color info
    loops$colorNames<-NA
    if (length(loopTrackColors)==1) {
      loops$colorNames<-loopTrackColors
    }
    else {
      loops$colorNames<-loopTrackColors[fileIndex]
    }

    allLoops<-rbind(allLoops, loops)


  }

  #specific color info
  if (is.null(specialLoopColors) == FALSE) {
    if (is.null(specialLoops) == TRUE) {
      message("error: must specify specialLoops if using specialLoopColors")
    }
    else {
      if (length(specialLoopColors)==1) {
        allLoops$colorNames[allLoops$name %in% specialLoops]<-specialLoopColors
      }
      else {
        for (specialLoopName in specialLoops) {
          loopColorIndex<-which(sapply(specialLoops,`==`, specialLoopName))
          loopColorName<-specialLoopColors[loopColorIndex]
          allLoops$colorNames[allLoops$name == specialLoopName]<-loopColorName
        }
      }
    }
  }

  allLoops$trackType<-"loops"
  allLoops$start<-(allLoops$start1+allLoops$end1)/2
  allLoops$end<-(allLoops$start2+allLoops$end2)/2

  allLoops$start<-allLoops$start/1000000
  allLoops$end<-allLoops$end/1000000

  if (loop_orientation == "above") {
    curvature=-0.5
    ymin<-0
    ymax<-1
  }
  else {
    curvature=0.5
    ymin<-(-1)
    ymax<-0

  }

  allLoops<-data.table::as.data.table(allLoops)
  allLoops$y_center<-0
  allLoops[, y_min := y_center + curvature]
  allLoops[, y_max := y_center]

  allLoops<-data.frame(allLoops)
  allLoops$start1<-allLoops$start1/1000000
  allLoops$end1<-allLoops$end1/1000000
  allLoops$start2<-allLoops$start2/1000000
  allLoops$end2<-allLoops$end2/1000000
  allLoops$rasterize<-rasterizePlot
  allLoops$alpha<-alpha
  allLoops$lineSize<-lineSize
  if (is.null(minScore)==FALSE) { #any loops with a score below the minimum will be more transparent
    allLoops$alpha[is.na(allLoops$score)==FALSE & allLoops$score < minScore]<-alpha/10
    allLoops$lineSize[is.na(allLoops$score)==FALSE & allLoops$score < minScore]<-lineSize/2

  }
  allLoops$displayNames<-as.character(allLoops$displayNames)
  allLoops$displayNames<-factor(allLoops$displayNames, levels=loopTrackNames)

  p1<-ggplot2::ggplot(allLoops)+ggplot2::geom_curve(ggplot2::aes(x=.data$start, y=0, xend=.data$end, yend=0, color=.data$colorNames, alpha=.data$alpha, linewidth=.data$lineSize),
                                  curvature=curvature)+ggplot2::scale_alpha_identity()+ggplot2::scale_linewidth_identity()+
    ggplot2::geom_segment(ggplot2::aes(x=.data$start1, xend=.data$end1, y=0, yend=0,
                     color=.data$colorNames, alpha=.data$alpha),linewidth=3, lineend = "butt", linejoin = "mitre")+
    ggplot2::geom_segment(ggplot2::aes(x=.data$start2, xend=.data$end2, y=0, yend=0,
                     color=.data$colorNames,alpha=.data$alpha),linewidth=3, lineend = "butt", linejoin = "mitre")+
    ggplot2::scale_color_identity()+ggplot2::facet_grid(displayNames ~., switch="y", scales="free_y", drop=FALSE)+
    ggplot2::scale_y_continuous(limits = c(ymin, ymax),expand = c(0, 0))+
    ggplot2::coord_cartesian(xlim=c(as.numeric(start)/1000000, as.numeric(end)/1000000))+
    ggplot2::theme_classic()+ggplot2::xlab("Location (Mb)")+ggplot2::theme(axis.title.y=ggplot2::element_blank(),
                                                axis.text.y=ggplot2::element_blank(),
                                                axis.ticks.y=ggplot2::element_blank(),
                                                axis.line.y=ggplot2::element_blank(),
                                                plot.title=ggplot2::element_text(hjust=0.5),
                                                legend.position = "none",strip.text = ggplot2::element_text(size = fontSize))
  if (rasterizePlot==TRUE) {
    rast_layers <- lapply(p1$layers, function(layer) {
      if (inherits(layer$geom, "GeomCurve") || inherits(layer$geom,"GeomSegment")) {
        ggrastr::rasterise(layer, dpi = 300)
      } else {
        layer
      }
    })

    p1$layers<-rast_layers
  }


  return(p1)




}






#' plot gene models, coverage tracks, peaks, and/or loops across region of interest
#'
#' @param genomicLoc character vector. A gene name (ex:"ZFX") OR a genomic locus of the form c(chromosome, start coordinate, end coordinate) (ex: c("chrX", 24040226, 24232664))
#' @param mart a mart obtained by utilizing the useMart function in biomaRt
#' @param ensembl_set if 'mart' is not provided, which Ensembl annotation you want to pull from. Default is 'hsapiens_gene_ensembl'. Most ensembl sets follow this format; for example if you're using mouse it would be 'mmusculus_gene_ensembl'
#' @param gene_symbol nomenclature system you are using for genes. Default is 'hgnc_symbol'but switches to 'external_gene_name' for any ensembl_set that is not 'hsapiens_gene_ensembl'. If you are using other nomenclature systems (ex: 'mgi_symbol' for mouse), you must specify. If using a custom annotation file, you likely need to use 'gene_id'.
#' @param includeNPC boolean, whether to include genomic features that are not protein coding. Default is FALSE.
#' @param custom_anno file path to a .gtf/.gff if using a custom annotation rather than ensembl
#' @param covFiles character vector containing a list of file paths to any coverage files (typically in bigWig or bedGraph format)
#' @param covTrackNames character vector containing the label of each coverage track. Should be the same length as covFiles and in the same order. By default it will label them as "Cov_1", "Cov_2", etc
#' @param covTrackColors character vector containing the color of each coverage track. If only one color is provided, all tracks will be that color. If different tracks must be different colors, specify the color for each track in order. Default is "black"
#' @param peakFiles character vector containing a list of file paths to any peak files (.bed format)
#' @param peakTrackNames character vector containing the label of each peak track. Should be the same length as 'peakFiles' and in the same order. By default it will label them as "Peaks_1", "Peaks_2", etc
#' @param peakTrackColors character vector containing color of each peak track. If only one color is provided, all tracks will be that color. If different tracks must be different colors, specify the color for each track in order. Default is "black"
#' @param loopFiles character vector containing a list of file paths to any contact/loop files (.bedpe format)
#' @param loopTrackNames character vector containing the label of each loop track. Should be the same length as 'contactFiles' and in the same order. By default it will label them as "Loops_1", "Loops_2", etc
#' @param loopTrackColors character vector containing color of each loop track. If only one color is provided, all tracks will be that color. If different tracks must be different colors, specify the color for each track in order. Default is 'black'
#' @param includeGenome boolean. Whether to include gene models in the final figure. Almost always set to TRUE. Default is TRUE.
#' @param includeTranscripts boolean. Whether to plot transcripts. Default is FALSE.
#' @param includeTxtNames boolean. Whether to label transcript models. If set to TRUE, requires includeGenome and include Transcripts to be set to TRUE. Default is FALSE.
#' @param transcript_list list, used to specify the names of the transcripts you want to include (if you only want a specific set of transcripts). Requires includeTranscripts to be set to TRUE. overrides any specified transcript filters
#' @param supportedTranscriptsOnly boolean, whether to only include well-supported transcripts as defined by the 'transcript_filters' parameter. Default is TRUE.
#' @param transcript_filters character vector of which transcript filters to clear when supportedTranscriptsOnly is set to TRUE. Common choices from ensembl are 'transcript_gencode_basic', 'transcript_appris', and 'transcript_tsl'. the 'appris_options' and 'tsl_options' parameters allow you to specify which APPRIS classification or transcript support levels, respectively, that you want to keep. Any other filter will simply keep transcripts that have a non-empty field for that filter. Default is 'transcript_gencode_basic'
#' @param appris_options character vector of APPRIS filters to implement if 'transcript_appris' is provided in 'transcript_filters'. can specify both principal (P1-P5) and alternative transcripts (A1-A5). ex: c("P1","P2","A1") will give you principal 1, principal2, and alternative 1 transcripts. If 'transcript_appris' is provided as 'transcript_filter' and 'appris_options' is not specified, all principal transcripts will be kept
#' @param tsl_options character vector of transcript support levels you want to keep if 'transcript_tsl' is provided in 'transcript_filters'
#' @param tag_options character vector of which tags to keep if using a custom annotation and setting supportedTranscriptsOnly to TRUE.
#' @param canonicalTranscriptOnly boolean, whether to only show the canoncial transcript. Default is FALSE.
#' @param upDown a character vector of length 2. if genomicLoc is a gene, how many base pairs upstream (first value) and downstream (second value) of the gene you want in your figure. Default is c(2000,2000)
#' @param genome_position if includeGenome is set to TRUE but 'genome' is not in 'trackOrder_type', where do you want gene models in your figure. set to "bottom" if you want the genome to be at the bottom of the figure; "top" if you want it to be at the top. Default is 'bottom'
#' @param ymin numeric, minimum value for y-axis for coverage plot
#' @param ymax numeric, maximum value for y-axis for coverage plot
#' @param fillArea boolean, set this to FALSE if you do not want to fill in the area under your coverage tracks. Default is TRUE
#' @param logScale boolean, set this to TRUE if you want data plotted on a log scale. Default is FALSE
#' @param rasterizeCovPlot boolean, set this to TRUE if you want to rasterize the coverage plot. This can be useful as coverage files can be rather large, and that can mess up the plot when you save it as a vector file. Default is FALSE
#' @param labelAllPeaks boolean, set to TRUE if each individual peak should be named. assumes names are provided in the fourth column of a given .bed file
#' @param specialPeaks character vector listing the names of specific peaks you want to either label or put in a different color
#' @param labelSpecialPeaks boolean, set to TRUE if each special peak should be named. requires 'specialPeaks' to be specified
#' @param specialPeakColors character vector containing the color(s) you want to give your special peaks. requires 'specialPeaks' to be specified
#' @param labelStrand boolean, set to TRUE if the strand of each peak should be indicated. assumes strand information is in the sixth column of a given .bed file
#' @param strandColors character vector of length 2. First is the color for the + strand, the second is the color for the - strand. overrides 'peakTrackColors' but does NOT override 'specialPeakColors'
#' @param lineSize numeric, width of lines used to draw loops. Default is 0.8
#' @param alpha numeric, alpha value (transparency) for loops. Default is 0.8
#' @param minScore numeric; any loops with scores lower than this will be thinner and more transparent, making the loops that pass this threshold more visible
#' @param rasterizeLoopPlot boolean, set this to TRUE if you want to rasterize the loop plot. This can be useful as loop files can be rather large, and that can mess up the plot when you save it as a vector file. Default is FALSE
#' @param specialLoops character vector listing the names of specific loops you want to put in a different color
#' @param specialLoopColors character vector containing the color(s) you want to give your special loops. requires 'specialLoops' to be specified
#' @param loop_orientation either "above" or "below". "above" will draw loops above a horizontal axis; "below" will draw loops below a horizontal axis. Default is "above'
#' @param trackOrder_type character vector with the desired order of track types. Default is c("coverage","peaks","loops", "genome")
#' @param fontSize a numeric for desired font size. Default is 9.
#' @param saveFigure boolean, whether you want to save the final figure to a file. Default is FALSE.
#' @param figureName if 'saveFigure' is TRUE, what you want the filename to be. Default is "plot_genomic_tracks_figure"
#' @param figureFormat if 'saveFigure' is TRUE, what format you want to save it as. Default is 'png'
#' @return list of length 6: figure (patchwork object with all tracks); coveragePlot (ggplot of coverage tracks); peakPlot (ggplot of peak tracks); loopPlot (ggplot of loop tracks); genePlot (ggplot of gene/transcript models); plotTitle (title of plot)
#' @export


plot_genomic_tracks<-function(genomicLoc, includeGenome=TRUE,includeTranscripts=FALSE,
                              includeTxtNames=TRUE, transcript_list=NULL,
                              supportedTranscriptsOnly=TRUE,
                              transcript_filters=c("transcript_gencode_basic"),
                              appris_options=NULL, tsl_options=NULL, tag_options=NULL,
                              canonicalTranscriptOnly=FALSE,
                              mart=NULL, ensembl_set="hsapiens_gene_ensembl",
                              gene_symbol="hgnc_symbol", includeNPC=FALSE,custom_anno=NULL,
                              upDown=c(2000,2000), genome_position="bottom",
                              covFiles=NULL, covTrackNames=NULL,covTrackColors="black",ymin=NULL, ymax=NULL, fillArea=TRUE, logScale=FALSE, rasterizeCovPlot=FALSE,
                              peakFiles=NULL, peakTrackNames=NULL,
                              labelAllPeaks=FALSE, peakTrackColors="black",
                              specialPeaks=NULL, labelSpecialPeaks=FALSE, specialPeakColors=NULL,
                              labelStrand=FALSE, strandColors=NULL,
                              loopFiles=NULL, loopTrackNames=NULL,
                              loopTrackColors="black", lineSize=0.8,alpha=0.8,minScore=NULL,rasterizeLoopPlot=FALSE,
                              specialLoops=NULL,specialLoopColors=NULL, loop_orientation="above",
                              trackOrder_type=c("coverage","peaks","loops","genome"), fontSize=9,
                              saveFigure=FALSE, figureName="plot_genomic_tracks_figure",figureFormat="png") {
  fontSize_text<-fontSize/ggplot2::.pt
  #basic checks
  if (length(genomicLoc) !=1 & length(genomicLoc)!=3) {
    stop("error: Please provide either a gene name or set of genomic coordinates")
  }
  for (tType in trackOrder_type) {
    if (tType %in% c("coverage","peaks","loops","genome")==FALSE) {
      stop("error: unrecognized track type in trackOrder_type")
    }
  }

  if (is.null(covTrackNames)==FALSE | is.null(peakTrackNames)==FALSE | is.null(loopTrackNames)==FALSE) {
    allTrackNames<-c(covTrackNames, peakTrackNames, loopTrackNames)
    if (length(allTrackNames) != length(unique(allTrackNames))) {
      print("warning: duplicate names detected between tracks. may cause issues if output is run through trackDJ")
    }
  }
  if (includeGenome==TRUE | length(genomicLoc) ==1 | "genome" %in% trackOrder_type) {
    a<-genomicLoc
    b<-includeTranscripts
    c<-transcript_list
    d<-supportedTranscriptsOnly
    e<-transcript_filters
    f<-appris_options
    g<-tsl_options
    h<-tag_options
    i<-canonicalTranscriptOnly
    j<-mart
    k<-ensembl_set
    l<-gene_symbol
    m<-includeNPC
    n<-custom_anno
    o<-upDown
    p<-includeTxtNames
    q<-fontSize

    gene_plot<-plot_gene(genomicLoc=a, includeTranscripts=b, includeTxtNames=p,transcript_list=c,
                         supportedTranscriptsOnly = d, transcript_filters = e,
                         appris_options =f,tsl_options=g, tag_options = h,
                         canonicalTranscriptOnly = i, mart=j, ensembl_set = k,
                         gene_symbol=l, includeNPC=m, custom_anno=n, upDown=o, fontSize=q)

    gene_plot<-gene_plot$genePlot
    gene_plot_data<-gene_plot$data
    gene_num<-base::nrow(gene_plot_data[is.na(gene_plot_data$type)==FALSE &  gene_plot_data$type=="Gene",])
    transcript_num<-base::nrow(gene_plot_data[is.na(gene_plot_data$type)==FALSE & gene_plot_data$type=="Transcript",])
    start<-unique(gene_plot_data$plot_start)
    end<-unique(gene_plot_data$plot_end)

    if ("chromosome_name" %in% colnames(gene_plot_data)) {
      chr<-unique(gene_plot_data$chromosome_name)

    }
    else {
      chr<-unique(gene_plot_data$chr)
    }

    chr<-as.character(chr)

    if (grepl("chr",chr)==TRUE) {
      chr<-substr(chr, 4, nchar(chr))
    }


    graphTitle<-get_plot_title(gene_plot)
    gene_plot<-gene_plot+ggplot2::ggtitle(NULL)
  }
  else {
    gene_plot<-NULL
    gene_num<-0
    transcript_num<-0
    if (length(genomicLoc)==3) { #if coordinates are already provided, we're golden
      chr<-as.character(genomicLoc[1])
      if (grepl("chr",chr)==TRUE) {
        chr<-substr(chr, 4, nchar(chr))
      }
      start<-as.numeric(genomicLoc[2])
      end<-as.numeric(genomicLoc[3])

      graphTitle<-paste0("chr",chr,":",start,"-",end)
    }
  }
  if ("genome" %in% trackOrder_type==FALSE & includeGenome==TRUE) {
    if (genome_position=="bottom" | genome_position=="below") {
      trackOrder_type<-c(trackOrder_type, "genome")
    }
    else {
      trackOrder_type<-c("genome",trackOrder_type)
    }
  }

  #figure out last one
  trackOrder_rev<-rev(trackOrder_type)
  lastType<-NULL
  for (typeName in trackOrder_rev) {
    if (typeName=="coverage") {
      checkTracks<-is.null(covFiles)
    }
    else if (typeName=="peaks") {
      checkTracks<-is.null(peakFiles)
    }
    else if (typeName=="loops" | typeName=="contacts") {
      checkTracks<-is.null(loopFiles)
    }
    else if (typeName=="genome") {
      if (includeGenome==TRUE) {
        checkTracks<-FALSE
      }
      else {
        checkTracks<-TRUE
      }
    }
    if (checkTracks == FALSE) {
      if (is.null(lastType)==TRUE) {
        lastType<-typeName
      }
    }
  }

  list_all<-list()
  size_all<-list()
  plotNum<-1

  for (tType in trackOrder_type) {
    if (tType=="coverage") {
      if (is.null(covFiles)==FALSE) {
        a<-c(chr,start,end)
        b<-covFiles
        c<-covTrackNames
        d<-covTrackColors
        e<-ymin
        f<-ymax
        g<-logScale
        h<-rasterizeCovPlot
        j<-fillArea
        k<-fontSize
        p_coverage<-plot_coverage(genomicLoc=a, covFiles=b, covTrackNames=c, covTrackColors=d, ymin=e, ymax=f, fillArea=j, logScale=g, rasterizePlot=h, fontSize=k)+ggplot2::theme(text = ggplot2::element_text(family = "Helvetica"))
        #if plot is not the last one, remove X axis
        if (lastType !="coverage") {
          p1<-p_coverage+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                               axis.ticks.x=ggplot2::element_blank(),
                               axis.line.x=ggplot2::element_blank(),
                               axis.text.x=ggplot2::element_blank())
        }
        else {
          p1<-p_coverage+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize),
                               axis.title.x=ggplot2::element_text(size=fontSize),
                               axis.line.x=ggplot2::element_line(linewidth=1))
        }
        p1<-p1+ggplot2::ggtitle(NULL)


        list_all[[plotNum]]<-p1
        plotNum<-plotNum+1

        newSize<-length(covFiles)
        size_all<-c(size_all,newSize)
      }
    }
    else if (tType=="peaks") {
      if (is.null(peakFiles)==FALSE) {
        a<-c(chr, start, end)
        b<-peakFiles
        c<-peakTrackNames
        d<-peakTrackColors
        e<-specialPeaks
        f<-labelSpecialPeaks
        g<-specialPeakColors
        h<-labelStrand
        i<-strandColors
        j<-labelAllPeaks
        k<-fontSize
        p_peaks<-plot_peaks(genomicLoc=a, peakFiles=b, peakTrackNames=c, peakTrackColors=d,
                            specialPeaks=e, labelSpecialPeaks=f, specialPeakColors=g,
                            labelStrand=h, strandColors=i, labelAllPeaks=j, fontSize=k)+ggplot2::theme(text = ggplot2::element_text(family = "Helvetica"))

        #if plot is not the last one, remove X axis
        if (lastType !="peaks") {
          p1<-p_peaks+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                            axis.ticks.x=ggplot2::element_blank(),
                            axis.line.x=ggplot2::element_blank(),
                            axis.text.x=ggplot2::element_blank())
        }

        else {
          p1<-p_peaks+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize),
                            axis.title.x=ggplot2::element_text(size=fontSize),
                            axis.line.x=ggplot2::element_line(linewidth=1))
        }

        p1<-p1+ggplot2::ggtitle(NULL)

        list_all[[plotNum]]<-p1
        plotNum<-plotNum+1

        newSize<-length(peakFiles)/2
        size_all<-c(size_all,newSize)
      }
    }

    else if (tType=="loops" |tType=="contacts") {
      if (is.null(loopFiles)==FALSE) {
        a<-c(chr, start, end)
        b<-loopFiles
        c<-loopTrackNames
        d<-loopTrackColors
        e<-lineSize
        f<-alpha
        g<-specialLoops
        h<-specialLoopColors
        i<-loop_orientation
        j<-minScore
        k<-rasterizeLoopPlot
        l<-fontSize

        p_loops<-plot_loops(genomicLoc=a, loopFiles=b, loopTrackNames=c, loopTrackColors=d, lineSize=e, alpha=f, specialLoops=g, specialLoopColors=h, loop_orientation=i, minScore = j, rasterizePlot=k, fontSize=l)+ggplot2::theme(text = ggplot2::element_text(family = "Helvetica"))

        #if plot is not the last one, remove X axis
        if (lastType !="loops" & lastType!="contacts") {
          p1<-p_loops+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                            axis.ticks.x=ggplot2::element_blank(),
                            axis.line.x=ggplot2::element_blank(),
                            axis.text.x=ggplot2::element_blank())
        }

        else {
          p1<-p_loops+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize),
                            axis.title.x=ggplot2::element_text(size=fontSize),
                            axis.line.x=ggplot2::element_line(linewidth=1))
        }

        p1<-p1+ggplot2::ggtitle(NULL)

        list_all[[plotNum]]<-p1
        plotNum<-plotNum+1

        newSize<-length(loopFiles)
        size_all<-c(size_all,newSize)

      }
    }

    else if (tType=="genome") {
      gene_plot<-gene_plot+ggplot2::theme(text = ggplot2::element_text(family = "Helvetica"))
      #if plot is not the last one, remove X axis
      if (lastType !="genome") {

        p1<-gene_plot+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                            axis.ticks.x=ggplot2::element_blank(),
                            axis.line.x=ggplot2::element_blank(),
                            axis.text.x=ggplot2::element_blank())

      }

      else {
        p1<-gene_plot+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize),
                            axis.title.x=ggplot2::element_text(size=fontSize),
                            axis.line.x=ggplot2::element_line(linewidth=1))
      }

      p1<-p1+ggplot2::ggtitle(NULL)
      list_all[[plotNum]]<-p1
      plotNum<-plotNum+1

      if (is.null(gene_plot$data)==FALSE) {
        geneSize<-length(unique(gene_plot$data$graph_location))/2
      }
      else {
        geneSize<-0
      }
      if (geneSize==0) {
        geneSize<-0.25
      }

      size_all<-c(size_all, geneSize)

    }

  }

  if (exists("p_coverage")==FALSE) {
    p_coverage<-NULL
  }
  if (exists("p_peaks")==FALSE) {
    p_peaks<-NULL
  }
  if (exists("p_loops")==FALSE) {
    p_loops<-NULL
  }

  #put them together

  finalFigure <- patchwork::wrap_plots(list_all,ncol = 1, heights = size_all)+
    patchwork:: plot_annotation(title = graphTitle,theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=fontSize)))

  if (saveFigure==TRUE) {
    if (substr(figureFormat,1,1)==".") {
      figureFormat<-substr(figureFormat,2,nchar(figureFormat))
    }
    size_total<-sum(size_all)

    ggplot2::ggsave(filename=paste0(figureName,".",figureFormat), finalFigure, device=figureFormat, dpi=320,bg="transparent", width=30, height=size_total*3, units="cm")
  }
  allPlots<-list(figure=finalFigure, coveragePlot=p_coverage, peakPlot=p_peaks, loopPlot=p_loops, genePlot=gene_plot, plotTitle=graphTitle)
  return(allPlots)

}

#' rearrange outputs from plot_genomic_tracks in the desired order
#'
#' @param plotList a list of outputs from plot_genomic_tracks
#' @param plotOrder character vector specifying the order of each track, using the names you provided to plot_genomic_tracks
#' @param fontSize a numeric for desired font size. Default is 9.
#' @param saveFigure boolean, whether you want to save the final figure to a file. Default is FALSE.
#' @param figureName if 'saveFigure' is TRUE, what you want the filename to be. Default is "trackDJ_figure"
#' @param figureFormat if 'saveFigure' is TRUE, what format you want to save it as. Default is 'png'
#' @return list of length two: figure (patchwork object with all tracks); singlePlots (list of ggplot objects of each individual track)
#' @export
trackDJ<-function(plotList, plotOrder, fontSize=9, saveFigure=FALSE, figureFormat="png", figureName="trackDJ_figure") {
  fontSize_text<-fontSize/ggplot2::.pt
  #see if plotOrder has "genome"
  if ("genome" %in% plotOrder==FALSE) {
    message("genome not found in plotOrder -- will not display genome")
  }
  #check plotList
  if (inherits(plotList,"ggplot")) {
    stop("error: plotList should be a list, not a ggplot object")
  }
  else if (length(intersect(names(plotList), c("figure", "coveragePlot","peakPlot","loopPlot","genePlot")))!=0) {
    plotList<-list(plotList)
  }

  covData<-data.frame()
  peakData<-data.frame()
  loopData<-data.frame()
  geneData<-data.frame()
  for (plotNum in 1:length(plotList)) {
    aPlot<-plotList[[plotNum]]
    aPlotNames<-names(aPlot)
    if (length(aPlotNames[aPlotNames %in% c("figure", "coveragePlot","peakPlot", "loopPlot", "genePlot","plotTitle")]) !=6) {
      stop("error: unrecognized plot type. must use output from plot_genomic_tracks")
    }
    if (plotNum == 1) { #get title and x coordinates
      graphTitle<-aPlot$plotTitle
      plotCoord<-ggplot2::ggplot_build(aPlot$figure)$layout$panel_params[[1]]$x.rang
      plot_start<-plotCoord[1]
      plot_end<-plotCoord[2]
    }
    if (is.null(plotList[[plotNum]]$coveragePlot)==FALSE) {
      y_range<-plotList[[plotNum]]$coveragePlot$coordinates$limits$y
      covData_raw<-plotList[[plotNum]]$coveragePlot$data
      covData_raw$ymin<-y_range[1]
      covData_raw$ymax<-y_range[2]

      covData<-rbind(covData, covData_raw)
    }
    if (is.null(plotList[[plotNum]]$peakPlot)==FALSE) {
      peakData<-rbind(peakData,plotList[[plotNum]]$peakPlot$data )
    }

    if (is.null(plotList[[plotNum]]$loopPlot)==FALSE) {
      loopData<-rbind(loopData, plotList[[plotNum]]$loopPlot$data)
    }

    if (is.null(plotList[[plotNum]]$genePlot)==FALSE) {
      newGeneData<-plotList[[plotNum]]$genePlot$data
      hasTxts<-TRUE
      if("geneName" %in% colnames(newGeneData)==FALSE) {
        newGeneData$geneName<-newGeneData$label
        newGeneData$transcriptName<-NA
        newGeneData$txt_label_loc<-NA
        newGeneData$includeTxtNames<-FALSE
        hasTxts=FALSE
      }

      if (base::nrow(geneData)==0) {
        geneData<-newGeneData
      }
      else {
        commonCols<-intersect(colnames(geneData), colnames(newGeneData))
        geneData<-geneData[,commonCols]
        newGeneData<-newGeneData[,commonCols]
        geneData<-rbind(geneData, newGeneData)
      }

    }



  }
  geneData<-geneData[,colnames(geneData) !="graph_location",]
  geneData <- geneData[!duplicated(geneData), ]

  #reorganize geneData to remove redundancy and correct graph locations
  geneData$start<-geneData$start*1000000
  geneData$end<-geneData$end*1000000
  if (base::nrow(geneData) != 0) {
    if (nrow(geneData[geneData$type=="Transcript",])==0) {
      final_geneData<-assign_y_levels_genes(geneData)

    }
    else {
      final_geneData<-assign_y_levels_transcripts(geneData)
    }

    final_geneData$start<-final_geneData$start/1000000
    final_geneData$end<-final_geneData$end/1000000
    #there shouldn't be any transcripts not affiliated with a gene, but just in case:
    notIncluded<-geneData[is.na(geneData$transcriptName)==FALSE & geneData$transcriptName %in% final_geneData$transcriptName==FALSE,]
    if (base::nrow(notIncluded)>0) {
      graphLoc<-min(final_geneData$graph_location)-3
      for (txtName in unique(notIncluded$transcriptName)) {
        singleTranscript<-notIncluded[notIncluded$transcriptName ==txtName & notIncluded$type=="Transcript",]
        if (base::nrow(singleTranscript)>1) {
          singleTranscript<-singleTranscript[1,]
        }
        singleTranscript$graph_location<-graphLoc
        singleTranscript$txt_label_loc<-graphLoc-0.75

        allExons<-geneData[is.na(geneData$transcriptName)==FALSE& geneData$transcriptName ==txtName & geneData$type=="Exon",]
        allExons<-allExons[!duplicated(allExons[,c("start","end")])]
        allExons$graph_location<-graphLoc
        allExons$txt_label_loc<-graphLoc-0.75

        final_geneData<-rbind(final_geneData, singleTranscript,allExons)

        graphLoc<-graphLoc-2
      }
    }
  }



  #make individual tracks:
  list_all<-list()
  list_single<-list()
  size_all<-list()
  size_total<-0
  plotNum<-1
  plotNames<-c()

  for (plotName in plotOrder) {
    if (plotName %in% covData$displayNames) {
      data_plot<-covData[covData$displayNames==plotName,]
      data_plot$displayNames<-factor(data_plot$displayNames, levels=plotName)

      if (unique(data_plot$logTransform)==FALSE){
        if (unique(data_plot$fillArea)==TRUE) {
          pCov<-ggplot2::ggplot(data_plot, ggplot2::aes(x = .data$start, y = .data$score, color=.data$colorNames, fill=.data$colorNames)) +
            ggplot2::geom_line() + ggplot2::geom_area()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
            ggplot2::facet_grid(displayNames ~ ., switch = "y", scales = "free_y")+ggplot2::xlab("Location (Mb)")+
            ggplot2::coord_cartesian(xlim = c(plot_start, plot_end), ylim = c(min(data_plot$ymin), max(data_plot$ymax)))+ggplot2::scale_y_continuous(n.breaks=3)+
            ggplot2::theme(axis.title.y=ggplot2::element_blank(),strip.text = ggplot2::element_text(size = fontSize),axis.text.y=ggplot2::element_text(size=fontSize),text = ggplot2::element_text(family = "Helvetica"))
        }
        else {
          pCov<-ggplot2::ggplot(data_plot, ggplot2::aes(x = .data$start, xend=.data$end,y = .data$score, color=.data$colorNames, fill=.data$colorNames)) +
            ggplot2::geom_segment()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
            ggplot2::facet_grid(displayNames  ~ ., switch = "y", scales = "free_y")+ggplot2::xlab("Location (Mb)")+
            ggplot2::coord_cartesian(xlim = c(plot_start, plot_end), ylim = c(min(data_plot$ymin), max(data_plot$ymax)))+ggplot2::scale_y_continuous(n.breaks=3)+
            ggplot2::theme(axis.title.y=ggplot2::element_blank(),strip.text = ggplot2::element_text(size = fontSize),axis.text.y=ggplot2::element_text(size=fontSize),text = ggplot2::element_text(family = "Helvetica"))
        }

      }

      else if (unique(data_plot$logTransform==TRUE)) {
        #plot with log scale
        if (unique(data_plot$fillArea)==TRUE) {
          pCov<-ggplot2::ggplot(data_plot, ggplot2::aes(x = .data$start, y = .data$score, color=.data$colorNames, fill=.data$colorNames)) +
            ggplot2::geom_line() + ggplot2::geom_area()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
            ggplot2::facet_grid(displayNames  ~ ., switch = "y", scales = "free_y")+ggplot2::xlab("Location (Mb)")+ggplot2::theme(strip.text = ggplot2::element_text(size = fontSize))+
            ggplot2::scale_y_continuous(trans='log2',n.breaks = 3)+
            ggplot2::coord_cartesian(xlim = c(plot_start, plot_end), ylim = c(min(data_plot$ymin), max(data_plot$ymax)))+ggplot2::theme(axis.title.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=fontSize),text = ggplot2::element_text(family = "Helvetica"))
        }
        else {
          pCov<-ggplot2::ggplot(data_plot, ggplot2::aes(x = .data$start, xend=.data$end, y = .data$score, color=.data$colorNames, fill=.data$colorNames)) +
            ggplot2::geom_segment()+ggplot2::scale_color_identity()+ggplot2::scale_fill_identity()+ggplot2::theme_classic()+
            ggplot2::facet_grid(displayNames  ~ ., switch = "y", scales = "free_y")+ggplot2::xlab("Location (Mb)")+ggplot2::theme(strip.text = ggplot2::element_text(size = fontSize))+
            ggplot2::scale_y_continuous(trans='log2', n.breaks=3)+
            ggplot2::coord_cartesian(xlim = c(plot_start, plot_end), ylim = c(min(data_plot$ymin), max(data_plot$ymax)))+ggplot2::theme(axis.title.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=fontSize),text = ggplot2::element_text(family = "Helvetica"))
        }

      }

      #rasterize coverage plot if necessary
      if (unique(data_plot$rasterize)==TRUE) {
        rast_layers <- lapply(pCov$layers, function(layer) {
          if (inherits(layer$geom, "GeomArea") || inherits(layer$geom, "GeomLine") || inherits(layer$geom, "GeomSegment")) {
            ggrastr::rasterise(layer, dpi = 300)
          } else {
            layer
          }
        })

        pCov$layers<-rast_layers
      }

      list_single[[plotNum]]<-pCov+ggplot2::ggtitle(graphTitle)+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line=ggplot2::element_line(linewidth=1),plot.title=ggplot2::element_text(hjust=0.5))
      #if plot is not the last one, remove X axis
      if (plotName != plotOrder[length(plotOrder)] | plotName %in% peakData$displayNames | plotName %in% loopData$displayNames) {
        pCov<-pCov+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank(),
                         axis.line.x=ggplot2::element_blank(),
                         axis.text.x=ggplot2::element_blank())

      }
      else {
        pCov<-pCov+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line.x=ggplot2::element_line(linewidth=1))
      }

      list_all[[plotNum]] <- pCov
      plotNum<-plotNum+1
      size_all<-c(size_all,1)
      size_total<-size_total+1
      plotNames<-c(plotNames,plotName)

    }

    if (plotName %in% peakData$displayNames) {
      data_plot<-peakData[peakData$displayNames==plotName,]
      data_plot$displayNames<-factor(data_plot$displayNames, levels=plotName)
      pPeak<-ggplot2::ggplot(data_plot, ggplot2::aes(x = .data$start, color=.data$colorNames)) +
        ggplot2::geom_segment(data=data_plot, ggplot2::aes(x=.data$start, xend=.data$end, y=1, yend=1,
                                         color=.data$colorNames),linewidth=3, lineend = "butt", linejoin = "mitre")+
        ggplot2::scale_color_identity()+ggplot2::theme_classic()+
        ggplot2::coord_cartesian(xlim = c(plot_start, plot_end), ylim=c(0,2))+
        ggplot2::facet_grid(displayNames ~ ., switch = "y", scales = "free_y")+
        ggplot2::xlab("Location (Mb)")+
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              axis.line.y = ggplot2::element_blank(),strip.text = ggplot2::element_text(size = fontSize),text = ggplot2::element_text(family = "Helvetica"))


      if (unique(data_plot$labelStrand) == TRUE) {
        pPeak<-pPeak+ggplot2::geom_text(data=data_plot, ggplot2::aes(x=.data$midX, label=.data$strand, y=.data$graph_location), size=fontSize_text, color="black", vjust = 0.5)

      }

      if (unique(data_plot$labelAllPeaks) == TRUE | unique(data_plot$labelSpecialPeaks)== TRUE) {
        pPeak<-pPeak+ggrepel::geom_text_repel(data=data_plot, y=0.85, ggplot2::aes(x=.data$midX, label=.data$label), size=fontSize_text, color="black", nudge_y=-0.2, direction="y", segment.size=0.5,na.rm = TRUE)
      }
      list_single[[plotNum]]<-pPeak+ggplot2::ggtitle(graphTitle)+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line.x=ggplot2::element_line(linewidth=1),plot.title=ggplot2::element_text(hjust=0.5))
      #if plot is not the last one, remove X axis
      if (plotName != plotOrder[length(plotOrder)]| plotName %in% loopData$displayNames) {
        pPeak<-pPeak+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank(),
                           axis.line.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank())
      }

      else {
        pPeak<-pPeak+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line=ggplot2::element_line(linewidth=1))
      }

      list_all[[plotNum]] <- pPeak
      plotNum=plotNum+1
      size_all<-c(size_all,0.5)
      size_total<-size_total+0.5
      plotNames<-c(plotNames,plotName)

    }

    if (plotName %in% loopData$displayNames) {
      data_plot<-loopData[loopData$displayNames==plotName,]
      data_plot$displayNames<-factor(data_plot$displayNames, levels=plotName)

      if (unique(data_plot$y_min) == (-0.5)) {
        curvature=-0.5
        ymin<-0
        ymax<-1
      }

      else {
        curvature=0.5
        ymin<-(-1)
        ymax<-0
      }

      pLoop<-ggplot2::ggplot(data_plot)+ggplot2::geom_curve(ggplot2::aes(x=.data$start, y=0, xend=.data$end, yend=0, color=.data$colorNames, alpha=.data$alpha, linewidth=.data$lineSize),
                                          curvature=curvature)+ggplot2::scale_alpha_identity()+ggplot2::scale_linewidth_identity()+
        ggplot2::geom_segment(ggplot2::aes(x=.data$start1, xend=.data$end1, y=0, yend=0,
                         color=.data$colorNames, alpha=.data$alpha),linewidth=3, lineend = "butt", linejoin = "mitre")+
        ggplot2::geom_segment(ggplot2::aes(x=.data$start2, xend=.data$end2, y=0, yend=0,
                         color=.data$colorNames,alpha=.data$alpha),linewidth=3, lineend = "butt", linejoin = "mitre")+
        ggplot2::scale_color_identity()+ggplot2::facet_grid(displayNames  ~., switch="y", scales="free_y")+
        ggplot2::scale_y_continuous(limits = c(ymin, ymax),expand = c(0, 0))+
        ggplot2::coord_cartesian(xlim=c(plot_start, plot_end))+
        ggplot2::theme_classic()+ggplot2::xlab("Location (Mb)")+ggplot2::theme(axis.title.y=ggplot2::element_blank(),
                                                    axis.text.y=ggplot2::element_blank(),
                                                    axis.ticks.y=ggplot2::element_blank(),
                                                    axis.line.y=ggplot2::element_blank(),
                                                    legend.position = "none",strip.text = ggplot2::element_text(size = fontSize),text = ggplot2::element_text(family = "Helvetica"))
      list_single[[plotNum]]<-pLoop+ggplot2::ggtitle(graphTitle)+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line.x=ggplot2::element_line(linewidth=1),plot.title=ggplot2::element_text(hjust=0.5))

      #if plot is not the last one, remove X axis
      if (plotName != plotOrder[length(plotOrder)]) {
        pLoop<-pLoop+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank(),
                           axis.line.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank())
      }
      else {
        pLoop<-pLoop+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line=ggplot2::element_line(linewidth=1))
      }
      #rasterize loop plot if necessary
      if (unique(data_plot$rasterize)==TRUE) {
        rast_layers <- lapply(pLoop$layers, function(layer) {
          if (inherits(layer$geom, "GeomCurve") || inherits(layer$geom, "GeomSegment")) {
            ggrastr::rasterise(layer, dpi = 300)
          } else {
            layer
          }
        })

        pLoop$layers<-rast_layers
      }


      list_all[[plotNum]] <- pLoop
      plotNum=plotNum+1
      size_all<-c(size_all,1)
      size_total<-size_total+1
      plotNames<-c(plotNames,plotName)

    }

    if (plotName=="genome") {
      arrow_properties_gene <- grid::arrow(type = "closed", angle = 20, length = grid::unit(0.3, "cm"))
      #if there are no transcripts:
      if (base::nrow(final_geneData[final_geneData$type=="Transcript",])==0) {
        pGene<-ggplot2::ggplot(final_geneData, ggplot2::aes(x = .data$start, y = .data$graph_location)) +
          ggplot2::geom_segment(ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                       data=final_geneData[final_geneData$type=="Gene",], linewidth=2.5, arrow=arrow_properties_gene, lineend = "butt",
                       linejoin = "mitre")+
          ggplot2::geom_label(color="black", fill="white",ggplot2::aes(x=.data$midX, y=.data$graph_location, label=.data$label),
                     fontface = "italic",size = fontSize_text)+
          ggplot2::xlim(plot_start, plot_end)+
          ggplot2::ylim(min(final_geneData$graph_location)-2,2)+ggplot2::theme_classic()+
          ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.line.y = ggplot2::element_blank(), plot.margin=ggplot2::margin(t=7.5, r=5.5,b=7.5, l=5.5, unit="pt"),text = ggplot2::element_text(family = "Helvetica"))+
          ggplot2::xlab("Location (Mb)")
      }

      #if there are transcripts:
      else {
        pGene<-ggplot2::ggplot(final_geneData, ggplot2::aes(x = .data$start, y = .data$graph_location)) +
          ggplot2::geom_segment(data=final_geneData[final_geneData$type=="Transcript",], ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                       linewidth=0.5, lineend = "butt", linejoin = "mitre")+
          ggplot2::geom_segment(data=final_geneData[final_geneData$type=="Exon",], ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                       linewidth=1.5, lineend = "butt", linejoin = "mitre")+
          ggplot2::geom_segment(data=final_geneData[final_geneData$type=="Gene",],ggplot2::aes(x=.data$start, xend=.data$end, y=.data$graph_location, yend=.data$graph_location),
                       linewidth=3, arrow=arrow_properties_gene, lineend = "butt", linejoin = "mitre")+
          ggplot2::geom_text(data=final_geneData[final_geneData$type=="Transcript" & final_geneData$includeTxtNames==TRUE,], color="black", ggplot2::aes(x=.data$midX,y=.data$txt_label_loc, label=.data$label), size=fontSize_text)+
          ggplot2::geom_label(data=final_geneData[final_geneData$type=="Gene",],color="black", fill="white",ggplot2::aes(x=.data$midX, y=.data$graph_location, label=.data$label),
                     fontface = "italic",size = fontSize_text)+
          ggplot2::xlim(plot_start, plot_end)+
          ggplot2::ylim(min(final_geneData$graph_location)-2,2)+ggplot2::theme_classic()+
          ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.line.y = ggplot2::element_blank(),text = ggplot2::element_text(family = "Helvetica"))+
          ggplot2::xlab("Location (Mb)")
      }
      list_single[[plotNum]]<-pGene+ggplot2::ggtitle(graphTitle)+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line.x=ggplot2::element_line(linewidth=1),plot.title=ggplot2::element_text(hjust=0.5))
      #if plot is not the last one, remove X axis
      if (plotName != plotOrder[length(plotOrder)]) {
        pGene<-pGene+ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank(),
                           axis.line.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank())
      }

      else {
        pGene<-pGene+ggplot2::theme(axis.text=ggplot2::element_text(size=fontSize), axis.title.x=ggplot2::element_text(size=fontSize), axis.line=ggplot2::element_line(linewidth=1))
      }

      gene_scale<-(base::nrow(final_geneData[final_geneData$type=="Gene",])/2)+(base::nrow(final_geneData[final_geneData$type=="Transcript",])/4)
      if(gene_scale==0) {
        gene_scale=0.25
      }
      list_all[[plotNum]] <- pGene
      plotNum<-plotNum+1
      size_all<-c(size_all,gene_scale)
      size_total<-size_total+gene_scale
      plotNames<-c(plotNames,plotName)


    }


  }

  #put them together
  p1 <- patchwork::wrap_plots(list_all,ncol = 1, heights = size_all)+  patchwork::plot_annotation(
    title = graphTitle,
    theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=fontSize))
  )
  if (length(unique(plotNames))!=length(plotNames)) {
    for (aName in unique(plotNames)) {
      if (length(plotNames[plotNames==aName])>1) {
        newNames<-paste0(aName,"_",1:length(plotNames[plotNames==aName]))
        plotNames[plotNames==aName]<-newNames
      }
    }
  }
  names(list_single)<-plotNames

  if (saveFigure==TRUE) {
    if (substr(figureFormat,1,1)==".") {
      figureFormat<-substr(figureFormat,2,nchar(figureFormat))
    }
    ggplot2::ggsave(filename=paste0(figureName,".",figureFormat), p1, device=figureFormat, dpi=320,bg="transparent", width=30, height=size_total*3, units="cm")
  }

  return(list(figure=p1, singlePlots=list_single))
}



#' obtain plot title from generated plot
#' @param p ggplot object
#' @return plot title

get_plot_title <- function(p) {
  if (!is.null(p$labels$title)) {
    return(p$labels$title)
  } else {
    return(NULL)   # no title set
  }
}

#' obtain list of chromosomes in a bigwig, bed, bedgraph, txt or tsv file
#' @param file filepath for bigwig, bed, bedgraph, txt or tsv file
#' @return chromosomes in file
get_chrom_names <- function(file) {


  ext <- tolower(tools::file_ext(file))

  #for bigwig files:
  if (ext %in% c("bw", "bigwig")) {
    si <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(file))
    return(GenomeInfoDb::seqnames(si))
  }

  # for bed and bedgraph files:
  if (ext %in% c("bed", "bedgraph", "bg", "txt", "tsv","bdg", "bedGraph")) {
    #only read first column (chromosome)
    chroms <- unique(data.table::fread(file, select = 1, header = FALSE)[[1]])
    return(chroms)
  }
  stop("Unsupported file type. Supported: bigwig (.bw), bed, bedGraph, bg, txt, tsv")
}




#' use import or tabular import for bedgraph files
#' @param file bedgraph file to read
#' @param chrom chromosome for the data you want to import
#' @param start_region start coordinate
#' @param end_region end coordinate
#' @return dataframe of data from bedgraph

read_bedgraph_region <- function(file, chrom, start_region, end_region) {
  bed <- NULL

  # Try importing with rtracklayer
  bed <- tryCatch({
    rtracklayer::import(file, format = "bedGraph")
  }, error = function(e) {
    message("Import failed: ", e$message)
    NULL
  })

  # If import failed, read manually
  if (is.null(bed)) {
    message("Reading manually as tab-delimited file.")
    bed_df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)

    if (ncol(bed_df) < 4) stop("BedGraph file must have at least 4 columns")
    colnames(bed_df)[1:4] <- c("chrom", "start", "end", "score")

    bed <- GenomicRanges::GRanges(
      seqnames = bed_df$chrom,
      ranges = IRanges::IRanges(start = bed_df$start + 1, end = bed_df$end),
      score = bed_df$score
    )
  }

  # Ensure score column exists
  if (is.null(S4Vectors::mcols(bed)$score)) {
    S4Vectors::mcols(bed)$score <- NA
  }

  # Define region of interest
  region <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = start_region, end = end_region))

  # Find overlaps
  hits <- GenomicRanges::findOverlaps(bed, region)
  bed_subset <- bed[S4Vectors::queryHits(hits)]

  # Convert to data frame
  df <- data.frame(
    chrom = as.character(GenomeInfoDb::seqnames(bed_subset)),
    start = start(bed_subset),
    end = end(bed_subset),
    score = S4Vectors::mcols(bed_subset)$score
  )

  return(df)
}

#' determine vertical coordinates of genes in gene model plot
#' @param df dataframe of gene chromosome, start coordinate, and end coordinate
#' @return dataframe complete with y coordinates for each gene in gene model plot

assign_y_levels_genes<-function(df) {
  chrName<-unique(df[,1])
  plottedSegments<-data.frame(start=min(df[1,]$start, df[1,]$end), end=max(df[1,]$start, df[1,]$end), loc=0)
  totalRows<-nrow(df)
  if (totalRows==1) {
    df$graph_location<-0
  }
  else {
    for (a in 2:totalRows) {
      geneStart<-min(df[a,]$start, df[a,]$end)
      geneEnd<-max(df[a,]$start, df[a,]$end)

      new_gr<-GenomicRanges::GRanges(seqnames=chrName, ranges=IRanges::IRanges(start=geneStart, end=geneEnd))

      if (a %% 2 ==0) { #evens can be on zero by default

        segList<-plottedSegments[plottedSegments$loc==0,]

        segList_gr<-GenomicRanges::GRanges(seqnames=chrName, ranges=IRanges::IRanges(start=segList$start, end=segList$end))
        hits<-GenomicRanges::findOverlaps(segList_gr, new_gr)

        if (length(hits)==0) {

          newLoc<-0
        }
        else {

          segList<-plottedSegments[plottedSegments$loc==(-3),]
          if (nrow(segList)==0) {
            newLoc<-(-3)
          }
          else {
            segList_gr<-GenomicRanges::GRanges(seqnames=chrName, ranges=IRanges::IRanges(start=segList$start, end=segList$end))
            hits<-GenomicRanges::findOverlaps(segList_gr, new_gr)
            if (length(hits)==0) {
              newLoc<-(-3)
            }
            else {
              b<-(-3)
              assigned<-FALSE
              while (assigned==FALSE) {
                b<-b-3
                if (b %in% plottedSegments$loc) {
                  segList<-plottedSegments[plottedSegments$loc==b,]
                  segList_gr<-GenomicRanges::GRanges(seqnames=chrName, ranges=IRanges::IRanges(start=segList$start, end=segList$end))

                  hits<-GenomicRanges::findOverlaps(segList_gr, new_gr)
                  if (length(hits)==0) {
                    newLoc<-b
                    assigned<-TRUE
                  }
                }
                else {
                  newLoc<-b
                  assigned<-TRUE
                }
              }
            }

          }


        }
      }
      else { #odds will be on negative three by default
        if ((-3) %in% plottedSegments$loc == FALSE) {
          newLoc<-(-3)
        }
        else {
          segList<-plottedSegments[plottedSegments$loc==(-3),]
          segList_gr<-GenomicRanges::GRanges(seqnames=chrName, ranges=IRanges::IRanges(start=segList$start, end=segList$end))
          hits<-GenomicRanges::findOverlaps(segList_gr, new_gr)
          if (length(hits)==0) {
            newLoc<-(-3)
          }
          else { #see if zero row works, if not make a new row
            segList<-plottedSegments[plottedSegments$loc==0,]
            segList_gr<-GenomicRanges::GRanges(seqnames=chrName, ranges=IRanges::IRanges(start=segList$start, end=segList$end))

            hits<-GenomicRanges::findOverlaps(segList_gr, new_gr)
            if (length(hits)==0) {
              newLoc<-0
            }
            else {
              b<-(-3)
              assigned<-FALSE
              while (assigned==FALSE) {
                b<-b-3
                if (b %in% plottedSegments$loc) {
                  segList<-plottedSegments[plottedSegments$loc==b,]
                  segList_gr<-GenomicRanges::GRanges(seqnames=chrName, ranges=IRanges::IRanges(start=segList$start, end=segList$end))

                  hits<-GenomicRanges::findOverlaps(segList_gr, new_gr)
                  if (length(hits)==0) {
                    newLoc<-b
                    assigned<-TRUE
                  }
                }
                else {
                  newLoc<-b
                  assigned<-TRUE
                }
              }
            }

          }
        }

      }
      newRow<-data.frame(start=geneStart, end=geneEnd, loc=newLoc)
      plottedSegments<-rbind(plottedSegments, newRow)
    }
    df$graph_location<-plottedSegments$loc
  }
  return(df)

}

#' determine vertical coordinates of transcripts in gene model plot
#' @param df dataframe of gene and transcript chromosome, start coordinate, and end coordinate
#' @return dataframe complete with y coordinates for each gene and transcript in gene model plot


assign_y_levels_transcripts<-function(df) {
  df<-df[,colnames(df) != "graph_location"]
  chrName<-unique(df[,1])
  gene_df<-df[df$type=="Gene",]

  geneLevels<-assign_y_levels_genes(gene_df)


  newDf<-geneLevels

  for (geneName in geneLevels$geneName) {
    loc<-unique(geneLevels[geneLevels$geneName==geneName,]$graph_location)

    txt_df<-df[df$geneName==geneName & df$type=="Transcript",]

    if (nrow(txt_df)>0) {
      for (a in 1:nrow(txt_df)) {
        loc<-loc-1.5
        singleTxt<-txt_df[a,]
        singleTxt$graph_location<-loc


        newDf<-rbind(newDf, singleTxt)

        txtName<-singleTxt$transcriptName


        exon_df<-df[df$type=="Exon" & df$transcriptName==txtName,]
        exon_df$graph_location<-loc
        newDf<-rbind(newDf, exon_df)
      }
    }

  }
  return(newDf)
}





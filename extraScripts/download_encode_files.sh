#!/bin/bash
# Download ENCODE files given their accession numbers.
# Usage: ./download_encode_files.sh ENCFF001ABC ENCFF002XYZ ...

# Check dependencies
command -v curl >/dev/null 2>&1 || { echo "❌ curl not found. Install it with: brew install curl"; exit 1; }
command -v jq >/dev/null 2>&1 || { echo "❌ jq not found. Install it with: brew install jq"; exit 1; }

BASE_URL="https://www.encodeproject.org"

for ACCESSION in "$@"; do
    echo "--------------------------------------"
    echo "🔍 Fetching metadata for: $ACCESSION"
    
    # Fetch JSON metadata
    JSON=$(curl -s -L -H "Accept: application/json" "${BASE_URL}/files/${ACCESSION}/")

    # Extract the download href and file name
    HREF=$(echo "$JSON" | jq -r '.href')
    FILENAME=$(basename "$HREF")

    if [[ "$HREF" == "null" || -z "$FILENAME" ]]; then
        echo "⚠️  Could not find file for accession $ACCESSION"
        continue
    fi

    DOWNLOAD_URL="${BASE_URL}${HREF}"
    echo "⬇️  Downloading ${FILENAME} ..."
    
    # Download the file
    curl -L -O "$DOWNLOAD_URL"

    echo "✅ Finished downloading $FILENAME"
done

echo "--------------------------------------"
echo "🎉 All downloads complete!"

#!/bin/bash

# Define input files
NL_RSEG_FILES=(
    "./RSEG/NIH-5943-CN-H3K9me3-WGL1-S1_S13_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
    "./RSEG/NIH-5533-CN-H3K9me3-WGL1-S1_S9_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
    "./RSEG/NIH-5577_CN_H3K9me3_S1_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
    "./RSEG/NIH-5444_CN_H3K9me3_S5_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
)

FXS_RSEG_FILES=(
    "./RSEG/NIH-5006_CN_H3K9me3_S13_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
    "./RSEG/NIH-5319_CN_H3K9me3_S9_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
    "./RSEG/NIH-5746_CN_H3K9me3_S1_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
    "./RSEG/NIH-6131_CN_H3K9me3_S13_ALTrBT2_SOMRD-DS_L-Filter-V2.bed"
)

BREACHES_FILE="./RSEG/BREACHes.bed"

# Combine all regions from the FXS set
cat "${FXS_RSEG_FILES[@]}" | sort -k1,1 -k2,2n | bedtools merge -i - > "$BREACHES_FILE"
echo "$BREACHES_FILE"

CONTROL_DOMAINS_FILE="./RSEG/control_domains.bed"
# Step 1: Combine all regions from the NL set of files
cat "${NL_RSEG_FILES[@]}" | sort -k1,1 -k2,2n | bedtools merge -i - -d 1000 | \
    awk '$3 - $2 > 25000' > "$CONTROL_DOMAINS_FILE"

# Step: Trim regions in BREACHES_FILE that overlap with CONTROL_DOMAINS_FILE
merged_BREACHES_FILE="./RSEG/BREACHes_merged.bed"
trimmed_BREACHES_FILE="./RSEG/BREACHes_trimmed.bed"

bedtools subtract -a <(sort -k1,1 -k2,2n "$BREACHES_FILE") \
                  -b <(sort -k1,1 -k2,2n "$CONTROL_DOMAINS_FILE") | \
    bedtools merge -i - -d 20000 > "$merged_BREACHES_FILE"

# Keep only regions larger than 1 Mb
awk '$3 - $2 > 1000000' "$merged_BREACHES_FILE" > "$trimmed_BREACHES_FILE"

rm "$BREACHES_FILE"

for ORIGINAL_FILE in "${FXS_RSEG_FILES[@]}"; do
    # Perform the overlap check: count the number of regions in the BREACHes_filtered file that overlap by >=10% with each of the original files
    OVERLAP_COUNT=$(bedtools intersect -a "$trimmed_BREACHES_FILE" -b "$ORIGINAL_FILE" -f 0.1 -u | wc -l)
    
    # Print the count for each original file to the console
    echo "Overlap count for $ORIGINAL_FILE: $OVERLAP_COUNT"
done

# Output file for the annotated BED file
ANNOTATED_BED_FILE="./annotated_BREACHES.bed"

# Create or clear the annotated BED file
> "$ANNOTATED_BED_FILE"

OUTPUT_TXT_FILE="upset_plot_data.txt"
> "$OUTPUT_TXT_FILE"

COLUMN_NAMES=()
for FILE in "${FXS_RSEG_FILES[@]}"; do
    BASENAME=$(basename "$FILE" .bed)
    COLUMN_NAMES+=("$BASENAME")
done

# Write column names as the header of the output file
echo -e "$(IFS=$'\t'; echo "${COLUMN_NAMES[*]}")" > "$OUTPUT_TXT_FILE"

# Initialize counters for each overlap count (1 to 4)
declare -A overlap_counts
for i in {1..4}; do
    overlap_counts[$i]=0
done

# Loop over each region BREACHes file
while read -r REGION; do
    OVERLAP_COUNT=0

    # Initialize a binary presence vector for upset plot data
    PRESENCE_VECTOR=()

    # Loop through each original BED file
    for ORIGINAL_FILE in "${FXS_RSEG_FILES[@]}"; do
        # Check if the current region overlaps with the original file by >=1%
        OVERLAP=$(echo "$REGION" | bedtools intersect -a <(echo "$REGION") -b "$ORIGINAL_FILE" -f 0.1 -u | wc -l)

        # Update presence vector and overlap count
        if [ "$OVERLAP" -gt 0 ]; then
            OVERLAP_COUNT=$((OVERLAP_COUNT + 1))
            PRESENCE_VECTOR+=(1)
        else
            PRESENCE_VECTOR+=(0)
        fi
    done

    echo -e "$(IFS=$'\t'; echo "${PRESENCE_VECTOR[*]}")" >> "$OUTPUT_TXT_FILE"

    # If the region overlaps with at least one file, append it with the overlap count to the new BED file
    if [ "$OVERLAP_COUNT" -gt 0 ]; then
        # Append the region along with the overlap count (annotated) to the new BED file
        echo -e "$REGION\t$OVERLAP_COUNT" >> "$ANNOTATED_BED_FILE"
        
        # Increment the corresponding overlap count
        if [ "$OVERLAP_COUNT" -le 4 ]; then
            overlap_counts[$OVERLAP_COUNT]=$((overlap_counts[$OVERLAP_COUNT] + 1))
        fi
    fi
done < "$trimmed_BREACHES_FILE"  # Read each region from the filtered BREACHes file

for ORIGINAL_FILE in "${FXS_RSEG_FILES[@]}"; do
    BASENAME=$(basename "$ORIGINAL_FILE" .bed)
    MERGED_FILE="./RSEG/${BASENAME}_BREACHes.bed"

    bedtools merge -i "$ORIGINAL_FILE" -d 10000 | \
    awk '$3 - $2 > 25000' | \
    bedtools intersect -b - -a "$trimmed_BREACHES_FILE" -f 0.1 -wa > "$MERGED_FILE"

done

rm "$CONTROL_DOMAINS_FILE"
rm "$trimmed_BREACHES_FILE"
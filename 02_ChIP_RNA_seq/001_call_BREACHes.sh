#!/bin/bash

# Define input files
FILTERED_BED_FILES_1=(
    "./RSEG/GMO7730_Domains_merged.bed"
    "./RSEG/GM09497_Domains_merged.bed"
    "./RSEG/FXS_Domains_merged.bed"
    "./RSEG/135_Domains_merged.bed"
)

ORIGINAL_BED_FILES=(
    "./RSEG/GM09497_Domains.bed"
    "./RSEG/FXS_Domains.bed"
    "./RSEG/GMO7730_Domains.bed"
    "./RSEG/135_Domains.bed"
)


# Merge regions and filter
for BED_FILE in ./RSEG/*Domains.bed; do
    # Get the filename without the directory or extension
    BASENAME=$(basename "$BED_FILE" .bed)
    MERGED_FILE="./RSEG/${BASENAME}_merged.bed"
    bedtools merge -i "$BED_FILE" -d 10000 | \
    awk '$3 - $2 > 250000' > "$MERGED_FILE"
done

BED_FILES_2=(
    "./RSEG/158_Domains_merged.bed"
    "./RSEG/176_Domains_merged.bed"
    "./RSEG/20b_Domains_merged.bed"
)

# Combine all regions from the first set of filtered files
COMBINED_FILE_1="./combined_all_regions_1.bed"
cat "${FILTERED_BED_FILES_1[@]}" | sort -k1,1 -k2,2n | bedtools merge -i - > "$COMBINED_FILE_1"

# Merge regions to create the BREACHes set (from the first set)
BREACHES_FILE="./BREACHes.bed"
bedtools merge -i "$COMBINED_FILE_1" -d 25000 > "$BREACHES_FILE"

rm "$COMBINED_FILE_1"

# Combine all regions from the second set of files 
COMBINED_FILE_2="./RSEG/combined_all_regions_2.bed"
cat "${BED_FILES_2[@]}" | sort -k1,1 -k2,2n | bedtools merge -i - > "$COMBINED_FILE_2"

CONTROL_DOMAINS_FILE="./RSEG/control_domains.bed"
bedtools merge -i "$COMBINED_FILE_2" -d 10000 > "$CONTROL_DOMAINS_FILE"

rm "$COMBINED_FILE_2"

#Remove regions in BREACHES that overlap with control_domains 
filtered_BREACHES_FILE="./BREACHes_filtered.bed"
bedtools subtract -a "$BREACHES_FILE" -b "$CONTROL_DOMAINS_FILE" | \
    bedtools merge -i - -d 100000 | \
    awk '$3 - $2 > 250000' > "$filtered_BREACHES_FILE"

rm "$BREACHES_FILE"

# Output file for the annotated BED file
ANNOTATED_BED_FILE="annotated_BREACHES.bed"

# Create or clear the annotated BED file
> "$ANNOTATED_BED_FILE"

OUTPUT_TXT_FILE="upset_plot_data.txt"

COLUMN_NAMES=()
for FILE in "${ORIGINAL_BED_FILES[@]}"; do
    BASENAME=$(basename "$FILE" .bed)
    COLUMN_NAMES+=("$BASENAME")
done

# Write column names as the header of the output file
echo -e "$(IFS=$'\t'; echo "${COLUMN_NAMES[*]}")" > "$OUTPUT_TXT_FILE"

# Initialize overlap counts for summary
declare -A overlap_counts
for i in {1..4}; do
    overlap_counts[$i]=0
done

# Loop over each region in the filtered BREACHes file
while read -r REGION; do
    OVERLAP_COUNT=0
    PRESENCE_VECTOR=()

    # Loop through each original BED file
    for ORIGINAL_FILE in "${FILTERED_BED_FILES_1[@]}"; do
        # Check if the current region overlaps with the original file by >=20%
        OVERLAP=$(echo "$REGION" | bedtools intersect -a <(echo "$REGION") -b "$ORIGINAL_FILE" -f 0.20 -u | wc -l)

        # Update presence vector and overlap count
        if [ "$OVERLAP" -gt 0 ]; then
            OVERLAP_COUNT=$((OVERLAP_COUNT + 1))
            PRESENCE_VECTOR+=(1)
        else
            PRESENCE_VECTOR+=(0)
        fi
    done

    # Append region with overlap count and presence vector to the annotated file
    echo -e "$REGION\t$OVERLAP_COUNT" >> "$ANNOTATED_BED_FILE"
    echo -e "$(IFS=$'\t'; echo "${PRESENCE_VECTOR[*]}")" >> "$OUTPUT_TXT_FILE"

    # Update summary counts
    if [ "$OVERLAP_COUNT" -le 5 ]; then
        overlap_counts[$OVERLAP_COUNT]=$((overlap_counts[$OVERLAP_COUNT] + 1))
    fi
done < "$filtered_BREACHES_FILE"

for ORIGINAL_FILE in "${FILTERED_BED_FILES_1[@]}"; do
    BASENAME=$(basename "$ORIGINAL_FILE" .bed)
    MERGED_FILE="./RSEG/${BASENAME}_BREACHes.bed"

    bedtools merge -i "$ORIGINAL_FILE" -d 20000 | \
    bedtools intersect -a "$filtered_BREACHES_FILE" -b - | \
    awk '$3 - $2 > 100000' > "$MERGED_FILE"

done

rm "$filtered_BREACHES_FILE"
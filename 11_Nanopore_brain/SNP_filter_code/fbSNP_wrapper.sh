mkdir -p filtered_bams

for bam in bam/*.bam; do
    python filter_by_SNP.py --bam "$bam" --out_dir filtered_bams
done

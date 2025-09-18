import pysam
from collections import Counter
import argparse
import os

def get_base_at_pos(read, ref_pos_0based):
    """Get base at 0-based reference position"""
    for query_pos, ref_aln_pos in read.get_aligned_pairs(matches_only=False):
        if ref_aln_pos == ref_pos_0based and query_pos is not None:
            return read.query_sequence[query_pos]
    return None

def compute_consensus(bam_path, chrom, positions_1based):
    """Compute consensus at 1-based positions"""
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # Convert to 0-based for internal processing
    positions_0based = [pos - 1 for pos in positions_1based]
    base_counts = {pos: Counter() for pos in positions_1based}
    
    # Fetch reads in the region (pysam.fetch uses 0-based coordinates)
    for read in bam.fetch(chrom, min(positions_0based), max(positions_0based) + 1):
        if read.is_unmapped or read.is_duplicate:
            continue
        
        for pos_1based, pos_0based in zip(positions_1based, positions_0based):
            base = get_base_at_pos(read, pos_0based)
            if base and base in "ACGT":
                base_counts[pos_1based][base] += 1
    
    bam.close()
    
    consensus_info = {}
    for pos in positions_1based:
        counts = base_counts[pos]
        total = sum(counts.values())
        consensus_base, count = counts.most_common(1)[0] if counts else ('N', 0)
        freqs = {base: count/total for base, count in counts.items()} if total > 0 else {}
        
        consensus_info[pos] = {
            "consensus": consensus_base,
            "counts": dict(counts),
            "freqs": freqs,
            "total": total,
        }
    
    return consensus_info

def filter_reads(bam_path, out_bam_path, chrom, consensus_bases_1based):
    """Filter reads based on consensus bases at 1-based positions"""
    bam = pysam.AlignmentFile(bam_path, "rb")
    out_bam = pysam.AlignmentFile(out_bam_path, "wb", header=bam.header)
    
    # Convert to 0-based for internal processing
    positions_1based = list(consensus_bases_1based.keys())
    positions_0based = [pos - 1 for pos in positions_1based]
    
    for read in bam.fetch(chrom, min(positions_0based), max(positions_0based) + 1):
        if read.is_unmapped or read.is_duplicate:
            continue
        
        match_all = True
        for pos_1based, pos_0based in zip(positions_1based, positions_0based):
            base = get_base_at_pos(read, pos_0based)
            if base != consensus_bases_1based[pos_1based]:
                match_all = False
                break
        
        if match_all:
            out_bam.write(read)
    
    bam.close()
    out_bam.close()
    pysam.index(out_bam_path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--out_dir", default="filtered_bams")
    parser.add_argument("--chrom", default="chrX")
    parser.add_argument("--positions", nargs="+", type=int, default=[147914300, 147917078])
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    sample = os.path.splitext(os.path.basename(args.bam))[0]
    out_bam_path = os.path.join(args.out_dir, f"{sample}.filtered.bam")
    out_summary_path = os.path.join(args.out_dir, f"{sample}.consensus_summary.txt")
    
    # Step 1: Compute consensus (positions are 1-based as input)
    consensus_info = compute_consensus(args.bam, args.chrom, args.positions)
    consensus_bases = {pos: consensus_info[pos]["consensus"] for pos in args.positions}
    
    # Step 2: Filter reads
    filter_reads(args.bam, out_bam_path, args.chrom, consensus_bases)
    
    # Step 3: Save summary
    with open(out_summary_path, "w") as f:
        for pos in args.positions:
            info = consensus_info[pos]
            f.write(f"Position: {args.chrom}:{pos}\n")
            f.write(f"  Total reads: {info['total']}\n")
            f.write(f"  Consensus base: {info['consensus']}\n")
            f.write(f"  Frequencies:\n")
            for base, freq in sorted(info["freqs"].items(), key=lambda x: -x[1]):
                f.write(f"    {base}: {freq:.3f}\n")
            f.write("\n")

if __name__ == "__main__":
    main()
# Example speedseq commands on a small slice of chromosome 20

# 1. Align with BWA
../bin/speedseq align \
    -o example \
    -M 3 \
    -p \
    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
    data/human_g1k_v37_20_42220611-42542245.fasta \
    data/NA12878.20slice.30X.fastq.gz

# 2. Detect SNVs and indels
../bin/speedseq var \
    -o example \
    data/human_g1k_v37_20_42220611-42542245.fasta \
    example.bam

# 3. Detect SVs
../bin/speedseq sv \
    -o example \
    -B example.bam \
    -S example.splitters.bam \
    -D example.discordants.bam \
    -R data/human_g1k_v37_20_42220611-42542245.fasta

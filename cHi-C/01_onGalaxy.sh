# cutadapt 1.16
cutadapt  -j ${GALAXY_SLOTS:-1}     -a 'TruSeq R1'='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'         -A 'TruSeq R2'='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'      --output='dataset1.gz' --paired-output='dataset2.gz'  --error-rate=0.1 --times=1 --overlap=3        --minimum-length=15 --pair-filter=any   --quality-cutoff=30   '2017_10_28_CB_CHIC_DL_S1_R_001_1.fq.gz' '2017_10_28_CB_CHIC_DL_S1_R_001_2.fq.gz'  > report.txt

# hicup 0.6.1
# bowtie2 2.2.6
# samtools 1.2
BOWTIE_PATH_BASH="$(which bowtie2)" && hicup_digester --re1 '^GATC' --genome 'mm10' mm10_UCSC.fa && mv *Digest_* digest_file.txt && hicup --zip --threads ${GALAXY_SLOTS:-1} --digest digest_file.txt --index 'mm10_UCSC' --bowtie2 $BOWTIE_PATH_BASH --keep    dataset1.gz dataset2.gz


# pysam 0.16.0 (for the DFL/PFL samples) or 0.13.0 (for the trunk samples)
# samtools 1.9
# This tool can be downloaded from https://testtoolshed.g2.bx.psu.edu/repository/download?repository_id=be5040251cd4afb7&changeset_revision=44365a4feb3b&file_type=gz
python fromHicupToJuicebox.py --fragmentFile digest_file.txt --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 --lineToSkipInFragmentFile 2 --methodForFrag hicup --useMid --output validPairs.txt hicup.bam && bash switchAndSort.sh validPairs.txt validPairs_sorted.txt

# Filtering with c10>=30 and c11>=30 with filtering tool from galaxy on validPairs_sorted.txt

# Filtering with (c3=='chr2' and c4<77000000 and c4>72402000) and (c7=="chr2" and c8<77000000 and c8>72402000) with filtering tool from galaxy on the previous output.txt

# cooler version 0.7.4
# tabix version 0.2.5
# h5py version 2.7.0

cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o validPairs_sorted_filtered.tabix validPairs_sorted_filtered.txt mm10.sizes

cooler makebins -o mm10.10kb.txt mm10.sizes 10000

cooler cload tabix --assembly mm10 -c2 7 -p2 8 mm10.10kb.txt validPairs_sorted_filtered.tabix output.cool

cooler balance --mad-max 5 --min-nnz 10 --min-count 0 --ignore-diags 2 --tol 1e-05 --max-iters 200 --cis-only -f output.cool

# Both tabix and cool files has been uploaded to GEO

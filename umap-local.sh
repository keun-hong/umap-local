#!/bin/sh

# Variables
in_folder=$1
in_genome=$2
kmer=$3

# Check if mandatory arguments are provided
if [ -z "$in_folder" ] || [ -z "$in_genome" ] || [ -z "$kmer" ]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 <input_folder_name> <genome_file_name> <kmer>"
    echo "Ex)    $0 data genome.fa '24 36 50 100'"
    echo ""
    exit 1
fi

date
echo ""
## ------------------------------------
# Variables
out_folder='mappability'
# Variables (auto-detected)
BOWTIEDIR=$(which bowtie|sed 's/\/bowtie//')
THREADS=$(expr $(nproc) - 2) # Max number of threads - 2

# Set short name pate
run_fold=$in_folder/$out_folder
# Find the minimum and maximum value
min=$(echo $kmer | tr ' ' '\n' | sort -n | head -1)
max=$(echo $kmer | tr ' ' '\n' | sort -n | tail -1)
# Genome file name
genome_name=$(echo $in_genome | cut -d '.' -f 1)
# Folders
WIGDIR=$run_fold/wgFiles
BEDDIR=$run_fold/bedFiles
BWDIR=$run_fold/bigWigFiles
LOGDIR=$run_fold/logFiles
mkdir -p $BEDDIR
mkdir -p $WIGDIR
mkdir -p $BWDIR
mkdir -p $LOGDIR
## ------------------------------------
# 0. Make chrom size file
echo "0. Make chrom size file"
in_genome_name=${in_genome%.fa}
in_chrsize=${in_genome_name}.chrom.sizes

if [ -f "${in_folder}/${in_genome}.gz" ]; then # decompress
    gzip -d "${in_folder}/${in_genome}.gz"
fi

samtools faidx ${in_folder}/${in_genome}
cut -f1,2 ${in_folder}/${in_genome}.fai > $in_folder/$in_chrsize
rm ${in_folder}/${in_genome}.fai

# 1. Run ubismap.py
echo "1. Run ubismap.py"
python ubismap.py $in_folder/$in_genome $in_folder/$in_chrsize $run_fold all.q $BOWTIEDIR/bowtie-build --kmer $kmer -write_script tmp.sh
rm tmp.sh
rm $in_folder/$in_chrsize

max_idx=$(expr $(tail -1 $run_fold/chrsize_index.tsv |awk '{ print $1 }') + 1) # Check max chr index

# 2. Generate bowtie index
echo "2. Generate bowtie index"

if [ -f "$run_fold/genome/Umap_bowtie.ind.1.ebwt" ]; then # Check if the bowtie index files already exist
    echo "Bowtie index already exists. Skipping index generation."
else
    bowtie-build $run_fold/genome/genome.fa $run_fold/genome/Umap_bowtie.ind \
        > $LOGDIR/02_index_genome.LOG 2> $LOGDIR/02_index_genome.ERR
fi

# 3. Bismap - GeNerate unique kmers
echo "3. Bismap - GeNerate unique kmers"
for k in $kmer;do
    echo "kmer: $k"
    parallel -j $THREADS python get_kmers.py $run_fold/chrsize.tsv $run_fold/kmers/k${k} $run_fold/chrs $run_fold/chrsize_index.tsv -job_id {1} --kmer k${k} \
        ">" $LOGDIR/03_Bismap.UniqueKmers_k${k}_{1}.LOG "2>" $LOGDIR/03_Bismap.UniqueKmers_k${k}_{1}.ERR ::: $(seq 1 $max_idx)
done

# 4. Bismap - Run Bowtie
echo "4. Bismap - Run Bowtie"
for k in $kmer; do
    echo "kmer: $k"
    parallel -j $THREADS python run_bowtie.py $run_fold/kmers/k${k} $BOWTIEDIR $run_fold/genome Umap_bowtie.ind -job_id {1} \
        ">" $LOGDIR/04_Bismap.RunBowtie_k${k}_{1}.LOG "2>" $LOGDIR/04_Bismap.RunBowtie_k${k}_{1}.ERR ::: $(seq 1 $max_idx)
done

# 5. Bismap - Unify bowtie outputs
echo "5. Bismap - Unify bowtie outputs"
for k in $kmer; do
    echo "kmer: $k"
    parallel -j $THREADS python unify_bowtie.py $run_fold/kmers/k${k} $run_fold/chrsize.tsv -job_id {1} \
        ">" $LOGDIR/05_Bismap.UnifyBowtie_k${k}_{1}.LOG "2>" $LOGDIR/05_Bismap.UnifyBowtie_k${k}_{1}.ERR ::: $(seq 1 $max_idx)
done

# 6. Move intermediary files
echo "6. Move intermediary files"
for k in $kmer; do
    echo "kmer: $k"
    mv $run_fold/kmers/k${k}/*kmer* $run_fold/kmers/k${k}/TEMPs \
        > $LOGDIR/06_Bismap.FileMov.LOG 2> $LOGDIR/06_Bismap.FileMov.ERR

    mv $run_fold/kmers/k${k}/*bowtie* $run_fold/kmers/k${k}/TEMPs \
        >> $LOGDIR/06_Bismap.FileMov.LOG 2>> $LOGDIR/06_Bismap.FileMov.ERR
done

# 7. Combine Umapped Files
echo "7. Combine Umapped Files"
for k in $kmer; do
    echo "kmer: $k"
    parallel -j $THREADS python combine_umaps.py $run_fold/kmers $run_fold/chrsize.tsv -job_id {1} \
        ">" $LOGDIR/07_Bismap.Combine_k${k}_{1}.LOG "2>" $LOGDIR/07_Bismap.Combine_k${k}_{1}.ERR ::: $(seq 1 $max_idx)
done

# 8. Generate BED files and wiggle files
echo "8. Generate BED files and wiggle files"
for k in $kmer; do
    echo "kmer: $k"
    parallel -j $THREADS python uint8_to_bed_parallel.py $run_fold/kmers/globalmap_k${min}tok${max} $run_fold/kmers/bedFiles ${genome_name}_umap -chrsize_path $run_fold/chrsize.tsv -bed -kmers $k -job_id {1} \
        ">" $LOGDIR/08_Bismap.uint8_to_bed_k${k}_{1}.LOG "2>" $LOGDIR/08_Bismap.uint8_to_bed_k${k}_{1}.ERR ::: $(seq 1 $max_idx)
    parallel -j $THREADS python uint8_to_bed_parallel.py $run_fold/kmers/globalmap_k${min}tok${max} $run_fold/kmers/wigFiles ${genome_name}_umap -chrsize_path $run_fold/chrsize.tsv -wiggle -kmers $k -job_id {1} \
        ">" $LOGDIR/08_Bismap.uint8_to_wig_k${k}_{1}.LOG "2>" $LOGDIR/08_Bismap.uint8_to_wig_k${k}_{1}.ERR ::: $(seq 1 $max_idx)
done

# 9. Merge the bed files and wiggle files of different chr
echo "9. Merge the bed files and wiggle files of different chr"
for k in $kmer; do
    echo "kmer: $k"
    python combine_wigs_or_beds.py $run_fold/kmers/bedFiles $BEDDIR --kmers $k \
        > $LOGDIR/09_Bismap.Combine_beds_k${k}_{1}.LOG 2> $LOGDIR/09_Bismap.Combine_beds_k${k}_{1}.ERR
    python combine_wigs_or_beds.py $run_fold/kmers/wigFiles $WIGDIR --kmers $k \
        > $LOGDIR/09_Bismap.Combine_wigs_k${k}_{1}.LOG 2> $LOGDIR/09_Bismap.Combine_wigs_k${k}_{1}.ERR
done

# 10. Convert the wiggle files to bigWig files
if [ ! -f ./wigToBigWig ]; then
    echo "wigToBigWig not found in the current directory. Downloading now..."
    # Download wigToBigWig from UCSC
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
    # Make the downloaded file executable
    chmod +x wigToBigWig
fi

files=$(ls $WIGDIR | grep wg)
for each in $files;do
    each2=$(echo $each | sed 's/wg.gz/wg/')
    newname=$(echo $each2 | sed 's/wg/bigWig/')
    gunzip $WIGDIR/$each
    ./wigToBigWig $WIGDIR/$each2 $run_fold/chrsize.tsv $BWDIR/$newname
    gzip -9 $WIGDIR/$each2
done

echo ""
echo "All tasks completed."
date
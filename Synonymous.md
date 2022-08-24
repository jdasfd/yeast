# Synonymous mutation among yeast

This markdown records the process of repeating the analysis from this [article](https://www.nature.com/articles/s41586-022-04823-w).

## Preparation

```bash
cd ~/Scripts
git clone https://github.com/wang-q/fig_table.git

cd fig_table
cp xlsx2csv.pl ~/data/yeast/scripts
```

- AmpUMI

```bash
cd ~/share
git clone https://github.com/pinellolab/AmpUMI.git
cd AmpUMI

vim AmpUMI.py
# add shebang line below:
#!/usr/bin/env python3
# save

python3 setup.py build
python3 setup.py install
sudo cp AmpUMI.py /usr/local/bin

cd /usr/local/bin
sudo chmod +x AmpUMI.py

AmpUMI.py -h
#usage: AmpUMI.py [-h] [--version] {Process,Collision,Distortion,CollisionNumber} ...
#
#AmpUMI - A toolkit for designing and analyzing amplicon sequencing experiments using unique molecular identifiers
```

- bbtools

```bash
brew install wang-q/tap/bbtools@37.77
```

## Usage intro

- `AmpUMI.py Process` is used for processing FASTQ reads after an amplicon sequencing experiment.

```bash
AmpUMI.py Process -h
```

```txt
usage: AmpUMI.py Process [-h] --fastq FASTQ --fastq_out FASTQ_OUT --umi_regex UMI_REGEX
                         [--min_umi_to_keep MIN_UMI_TO_KEEP] [--write_UMI_counts] [--write_alleles_with_multiple_UMIs]
```

- `bbmerge.sh` - merges overlapping or nonoverlapping pairs into a single reads.

```bash
bbmerge.sh --help
```

```txt
Description: Merges paired reads into single reads by overlap detection.
With sufficient coverage, can also merge nonoverlapping reads by kmer extension.
Kmer modes requires much more memory, and should be used with the bbmerge-auto.sh script.
Please read bbmap/docs/guides/BBMergeGuide.txt for more information.

Usage for interleaved files:    bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>
Usage for paired files:         bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>
```

## Seq files

All data were downloaded from the Bioproject [PRJNA750109](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA750109).

```bash
mkdir -p ~/data/yeast/ena
cd ~/data/yeast/ena

# download SraRunTable.txt to this doc
# PRJNA750109, search it in NCBI SRA Run Selector

# SraRunTable.txt is in .csv format
cat SraRunTable.txt |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    sed '1 s/^/#/' |
    keep-header -- tsv-sort -k2,2 -k3,3nr |
    tsv-uniq -H -f "Sample\ Name" --max 1 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

aria2c -j 4 -x 2 -s 2 --file-allocation=none -c -i ena_info.ftp.txt
# aria2c could be used, although it was better using aspera (shell need to be modified if using aspera)

md5sum --check ena_info.md5.txt
# check if there was any mistake during downloading
```

## Grouping seq files

### Get info

`mkdir -p ~/data/yeast/info`

Download all supplementary files and tables from the [website](https://www.nature.com/articles/s41586-022-04823-w) and save them into info dir.

```bash
mkdir -p ~/data/yeast/trim
cd ~/data/yeast

# get all DNA and RNA seq files
cat ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srx,srr |
    tsv-filter -H --regex 'name:\.t\d' |
    sed '1d' | 
    tr '.' '\t' |
    sed 's/t0/t0\t0/' |
    sed '1igene\ttime\tgroup\tsrx\tsrr' \
    > info/gene_time.tsv

# 41586_2022_4823_MOESM7_ESM.xlsx contains multiple sheets
perl scripts/xlsx2csv.pl -f info/41586_2022_4823_MOESM7_ESM.xlsx \
    --sheet "YPD DFE Sequencing primers " |
    sed '1,2d' |
    tr " " "_" |
    sed 's/__/_/' |
    sed 's/^"//' |
    sed 's/",/,/' | 
    perl -nla -F"," -e '
    $F[0] =~ /^(.+)_\d$/;
    $gene = $1;
    $F[1] =~ /(N+.+)$/;
    $primer = $1;
    print "$gene\t$primer";
    ' |
    uniq |
    sed 's/Reverse_i7/R/'\
    > info/DFE_seq_primers.tsv

cat info/DFE_seq_primers.tsv |
    cut -f 1 |
    cut -d '_' -f 1 |
    sort |
    uniq \
    > info/gene.lst
```

### Comparing SE and PE

According to the article, Illumina

```bash
mkdir ~/data/yeast/trim
cd ~/data/yeast/ena

# AmpUMI.py do not accept gz format
AmpUMI.py Process --fastq SRR15274411_1.fastq \
    --fastq_out ../trim/SRR15274411_1.dedup.fastq \
    --umi_regex "^NNNNNNNNAGACTTTAGGGCTCGGTAATT" \
    --write_UMI_counts \
    --write_alleles_with_multiple_UMIs

faops filter -l 0 SRR15273966_1.dedup.fastq stdout | grep '^>' | wc -l
#61461

faops filter -l 0 SRR15273966_2.dedup.fastq stdout | grep '^>' | wc -l
#16384

faops filter -l 0 SRR15274411_1.dedup.fastq stdout | grep '^>' | wc -l
#64655

faops filter -l 0 SRR15274411_2.dedup.fastq stdout | grep '^>' | wc -l
#16384
```

So for pair-end UMIs, two different files were slightly different after filtering by UMIs.

### Dealing with UMIs

First, using SE files to estimate fitness. How to deal with PE UMI-seq files will be studied afterwards.

- Grouping files

```bash
cd ~/data/yeast

# grouping and decompress for AmpUMI using
cat info/gene_time.tsv |
    sed '1d' |
    parallel --col-sep "\t" -k -j 4 '
    if ! [[ -d trim/repeat{3} ]]; then
        mkdir trim/repeat{3}
    fi
    cp ena/{5}_1.fastq.gz trim/repeat{3}/{1}_{2}_{3}.fastq.gz
    gzip -d trim/repeat{3}/{1}_{2}_{3}.fastq.gz
    '
```

- AmpUMI

Because `AmpUMI` only accepts decompressed files and does not accept input from stdin, so the input files are renamed and provided after grouping them.

The UMI sequences could be accepted by `AmpUMI`, so all the primers were included into the process.

```bash
cd ~/data/yeast

for gene in $(cat info/gene.lst)
do
    umi=$(cat info/DFE_seq_primers.tsv |
          tsv-filter --str-eq 1:${gene}_F |
          tsv-select -f 2)
    AmpUMI.py Process --fastq trim/repeat0/${gene}_t0_0.fastq \
        --fastq_out trim/repeat0/${gene}_t0_0.dedup.fastq \
        --umi_regex $umi \
        --write_UMI_counts \
        --write_alleles_with_multiple_UMIs
done
```

### Counting mutations

```bash

```

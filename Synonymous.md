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

## AmpUMI usage

`AmpUMI.py Process` is used for processing FASTQ reads after an amplicon sequencing experiment.

```bash
AmpUMI.py Process -h
```

```txt
usage: AmpUMI.py Process [-h] --fastq FASTQ --fastq_out FASTQ_OUT --umi_regex UMI_REGEX
                         [--min_umi_to_keep MIN_UMI_TO_KEEP] [--write_UMI_counts] [--write_alleles_with_multiple_UMIs]
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

```bash
mkdir -p /mnt/e/data/yeast/info
cd /mnt/e/data/yeast

# get all DNA and RNA seq files
cat ../ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srx,srr |
    tsv-filter -H --regex 'name:\.t\d' |
    sed 's/\./_/g' \
    > gene_time.tsv

xlsx2csv seq_primer.xlsx |
    sed 's/\s/_/g' |
    mlr --icsv --otsv cat \
    > seq_primer.tsv

cat seq_primer.tsv |
    cut -f 1 |
    cut -d '_' -f 1 |
    sort | uniq \
    > gene.lst

cat gene.lst |
    parallel -j 1 -k '
    echo "==> Trim {}"
    Fp=$(cat seq_primer.tsv | tsv-filter --str-eq 1:${}_amp_F | cut -f 2)
    Rp=$(cat seq_primer.tsv | tsv-filter --str-eq 1:${}_amp_R | cut -f 2)
    
    for file=$(cat gene_time.tsv | grep {} | cut -f 3)
    do
        trim_galore ../ena/${file}_1
    
    '
```

## Trimming

According to the article, Illumina 

```bash
mkdir ~/data/yeast/trim
cd ~/data/yeast/ena


```
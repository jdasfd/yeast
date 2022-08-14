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

## Info

`mkdir -p ~/data/yeast/info`

Download all supplementary files and tables from the [website](https://www.nature.com/articles/s41586-022-04823-w) and save them into info dir.

```bash
mkdir -p ~/data/yeast/info
cd ~/data/yeast

# get all DNA and RNA seq files
cat ../ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srx,srr |
    tsv-filter -H --regex 'name:\.t\d' |
    sed 's/\./_/g' \
    > gene_time.tsv

# 41586_2022_4823_MOESM7_ESM.xlsx contained multiple sheets
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
    uniq \
    > info/DFE_seq_primers.tsv

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

```bash
bbmerge.sh in1=SRR15273966_1.fastq in2=SRR15273966_2.fastq \
    interleaved=false out=SRR15273966.merge.fastq

faops size SRR15273966.merge.fastq |
    cut -f 2 |
    sort |
    uniq |
    tsv-summarize --max 1 --min 1
#309  101
```

It is meaning that the merge will not be automatically adapted to amplicon sequencing, so the method should be changed.

```bash
for gene in $(ls | perl -p -e 's/_.+$//')
do
    echo $gene
    cat ${gene}_genetic_interactions.txt |
        grep -v '^!' |
        grep -v '^\s*$' |
        sed '1d' |
        tsv-select -f 1,3,6,8 |
        tsv-filter --not-iregex 3:"Dosage" |
        tsv-filter --not-iregex 3:"Rescue" |
        tsv-summarize --group-by 4 --count
    echo
done
```

## Variants info

### Data preparation

- Extract vcf info from 1002 genomes project

```bash
mkdir ~/data/yeast/isolates
cd ~/data/yeast/isolates

wget http://1002genomes.u-strasbg.fr/isolates/page8/files/1002genomes.txt

# Remove references part
cat 1002genomes.txt |
    sed '1013,1042d' \
    > 1002genomes.tsv

# Human and Human, clinical were removed
for catgry in $(cat 1002genomes.tsv | sed '1d' | cut -f 4 | grep -v "^Human" | sort | uniq)
do
    echo "==> ${catgry}"
    cat 1002genomes.tsv |
        sed '1d' |
        tsv-select -f 1,2,4 |
        tsv-filter --iregex 3:${catgry} |
        tsv-select -f 1,2 \
        > ${catgry}.tsv
    cat ${catgry}.tsv |
        tsv-select -f 1 \
        > ${catgry}.lst
done
```

```bash
cd ~/data/yeast
mkdir vcf

# split vcf according to groups
for group in $(cat isolates/1002genomes.tsv | sed '1d' | cut -f 4 | grep -v "^Human" | sort | uniq)
do
    echo "==> ${group}"
    
    sample=$(cat isolates/${group}.lst |
        tr '\n' ',' |
        sed 's/,$//')
    
    bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz \
        --threads 8 -s ${sample} |
        bcftools +fill-tags -Ob -o vcf/1011Matrix.${group}.bcf
done
# -Ob: output bcf (compressed)
# only compressed bcf format could be index

cd ~/data/yeast/vcf

parallel -j 4 " \
bcftools index --threads 3 {} \
" ::: $(ls *.bcf)
```

- Extract all mutations in vcf-like format

```bash
cd ~/data/yeast
mkdir gene

perl scripts/xlsx2csv.pl -f info/41586_2022_4823_MOESM3_ESM.xlsx |
    sed '1d' |
    cut -d, -f 1 \
    > gene/stdname.lst

rm gene/std_sysname.tsv
for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"
    cat ../mrna-structure/sgd/saccharomyces_cerevisiae.gff |
        grep "${gene}" |
        cut -f 9 | 
        perl -nla -F';' -e '
        $F[1] =~ /Name=(.*)/ or next;
        $name = $1;
        $F[2] =~ /gene=(.*)/ or next;
        $gene = $1;
        print "$name\t$gene";
        ' \
        >> gene/std_sysname.tsv
done

wc -l gene/std_sysname.tsv
#34 gene/std_sysname.tsv
# more than 21 genes, need filtering

cat gene/std_sysname.tsv |
    tsv-select -f 2,1 |
    tsv-join --filter-file gene/stdname.lst \
    --key-fields 1 \
    > tmp.tsv && mv tmp.tsv gene/std_sysname.tsv

wc -l gene/std_sysname.tsv
#21 gene/std_sysname.tsv
# alright

cat gene/std_sysname.tsv |
    tsv-select -f 2 \
    > gene/sysname.lst

cat gene/std_sysname.tsv |
    tsv-select -f 2,1 \
    > gene/sys_stdname.tsv
# change standard names and system names
```

Locating mutations on the genome scale.

```bash
cd ~/data/yeast

perl scripts/xlsx2csv.pl -f info/41586_2022_4823_MOESM9_ESM.xlsx \
    --sheet "Fig. 2abc" |
    sed '1d' |
    tsv-select -d, -f 1,3,2 |
    mlr --icsv --otsv cat |
    tr '-' '\t' \
    > info/fit.tmp.tsv

for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"
    mkdir gene/${gene}
    cat info/fit.tmp.tsv |
        grep "^$gene" \
        > gene/${gene}/${gene}.tmp.tsv
done

# extract 150 coding regions (detected)
for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"
    cat gene/${gene}/${gene}.tmp.tsv |
        tsv-select -f 1,2,3 |
        tsv-uniq |
        sort -nk 2,2 |
        wc -l
    echo ">${gene}" \
        > gene/${gene}/${gene}.mut.fa
    cat gene/${gene}/${gene}.tmp.tsv |
        tsv-select -f 1,2,3 |
        tsv-uniq |
        sort -nk 2,2 | 
        perl -nae '
        chomp;
        $i = 1 if ($F[1] == 1);
        if($F[1] eq $i){
            print "$F[2]";
            $i++;
        }else{
            print "-";
            $i++;
            redo;
        }
        END{print "\n";}
        ' >> gene/${gene}/${gene}.mut.fa
done

# extract CDS region
# because of the coding region was selected from the SGD coding region
faops some ../mrna-structure/sgd/orf_genomic_all.fasta \
    gene/sysname.lst gene/21orf.fa

# rename by standard name
faops replace -l 0 gene/21orf.fa gene/sys_stdname.tsv gene/std_orf.fa

# align to the gene
for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"
    faops one -l 0 gene/std_orf.fa ${gene} gene/${gene}/${gene}.fa
    cat gene/${gene}/${gene}.fa gene/${gene}/${gene}.mut.fa |
        muscle -out gene/${gene}/${gene}.aln.fa -quiet
done

# extract CDS genomic location
cat ../mrna-structure/sgd/saccharomyces_cerevisiae.gff |
    perl -nla -e '
    next if /^#/;
    next if /^[ATCG]/;
    print "$F[0]\t$F[3]\t$F[4]\t$F[8]" if ($F[2] eq CDS);
    ' > gene/CDS.gff

for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"
    sys=$(cat gene/std_sysname.tsv |
              tsv-filter --str-eq 1:${gene} |
              tsv-select -f 2)
    cat gene/CDS.gff |
        tsv-filter --iregex 4:${sys} |
        perl -nla -e '
        $F[0] =~ /^chr(.+)/;
        $chr = $1;
        $F[3] =~ /Name=(.+)_CDS/;
        $name = $1;
        print "$name\t$chr\t$F[1]\t$F[2]";
        ' \
        > gene/${gene}/${gene}_region.tsv
done

# CDS location on the genomic region
for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"
    perl scripts/loc2vcf.pl \
        -r gene/${gene}/${gene}_region.tsv \
        -l gene/${gene}/${gene}.tmp.lst \
        -a gene/${gene}/${gene}.aln.fa \
        > gene/${gene}/${gene}.mut.vcf
done
```

- Filtering with bcf of different groups using region

```bash
mkdir ~/data/yeast/vcf/region
cd ~/data/yeast/vcf

echo -e "#CHROM\tPOS\tEND" > region/region.bed

for gene in $(cat ../gene/stdname.lst)
do
    cat ../gene/${gene}/${gene}_region.tsv |
        perl -nla -F"\t" -e '
            BEGIN {
                our %roman = (
                    "XVI"   => 16,
                    "XV"    => 15,
                    "XIV"   => 14,
                    "XIII"  => 13,
                    "XII"   => 12,
                    "XI"    => 11,
                    "X"     => 10,
                    "IX"    => 9,
                    "VIII"  => 8,
                    "VII"   => 7,
                    "VI"    => 6,
                    "V"     => 5,
                    "IV"    => 4,
                    "III"   => 3,
                    "II"    => 2,
                    "I"     => 1
                );
            }
            my $chr = $roman{$F[1]};
            print "chromosome$chr\t$F[2]\t$F[3]";
        ' \
    >> region/region.bed 
done

# transfer all mutation to 1 vcf
for gene in $(cat ../gene/stdname.lst)
do
    cat ../gene/${gene}/${gene}.mut.vcf \
        >> ../vcf/region/all.mut.vcf
done

cat region/all.mut.vcf | wc -l
#8548
# all mutation with fitness scores

for group in $(ls *.bcf | sed 's/^1011Matrix\.//' | sed 's/\.bcf$//')
do
bcftools view 1011Matrix.${group}.bcf -Ov -R region/region.bed \
    -o ${group}.vcf
done

for group in $(ls *.vcf | sed 's/\.vcf//')
do
    echo "==> ${group}"
    
    cat ${group}.vcf |
        perl -nla -F"\t" -e '
        /^\#\#/ and next;
        splice @F, 8;
        print join qq{\t}, @F;
    ' \
    > ${group}.tsv

    cat ${group}.tsv |
        perl -nla -F"\t" -e '
            BEGIN {
                our %roman = (
                    16 => "XVI",
                    15 => "XV",
                    14 => "XIV",
                    13 => "XIII",
                    12 => "XII",
                    11 => "XI",
                    10 => "X",
                    9  => "IX",
                    8  => "VIII",
                    7  => "VII",
                    6  => "VI",
                    5  => "V",
                    4  => "IV",
                    3  => "III",
                    2  => "II",
                    1  => "I"
                );
            }
            next if /^#/;
            my $loca = $F[0];
            $loca =~ /^chromosome(\d+)/;
            $chr = $roman{$1};
            my $R        = length $F[3];
            my $A        = length $F[4];
            my @info     = split /;/, $F[7];
            my @AF       = split /=/, $info[1];
            my $Freq_vcf = $AF[1];
            my @AC       = split /=/, $info[0];
            my @AN       = split /=/, $info[2];
            my $ALT_vcf  = $AC[1];
            my $REF_vcf  = $AN[1] - $AC[1];
    
            if ( $R == 1 && $A == 1 ) {
                print qq{$chr\t$F[1]\t$F[3]\t$F[4]\t$Freq_vcf\t$REF_vcf\t$ALT_vcf};
            }
        ' \
        > region/${group}.tsv
done
```

### Statistical analysis

- Join analysis

```bash
cd ~/data/yeast/vcf/region

for group in $(ls *.tsv | sed 's/\.tsv//')
do
    echo "==> ${group}"
    cat ${group}.tsv |
        tsv-filter --ne 5:0 |
        wc -l
    
    tsv-join all.mut.vcf -k 1,2,3 \
        --filter-file <(cat ${group}.tsv | tsv-filter --ne 5:0) |
        wc -l
done
```

==> Bakery
220
15
==> Beer
368
37
==> Bioethanol
152
18
==> Cider
135
12
==> Dairy
137
18
==> Distillery
246
32
==> Fermentation
237
27
==> Flower
138
12
==> Fruit
378
30
==> Industrial
233
18
==> Insect
249
15
==> Lab
80
12
==> Nature
535
42
==> Palm
274
27
==> Probiotic
58
0
==> Sake
143
18
==> Soil
328
27
==> Tree
480
37
==> Unknown
264
30
==> Water
214
12
==> strain
80
12
==> wine
596
60

- Extract fitness data

```bash
perl scripts/xlsx2csv.pl -f info/41586_2022_4823_MOESM9_ESM.xlsx \
    --sheet "Fig. 2abc" |
    sed '1d' |
    tsv-select -d, -f 1,3,2 |
    mlr --icsv --otsv cat |
    tr '-' '\t' \
    > info/fit.tmp.tsv

cd ~/data/yeast/vcf/region

for group in $(ls *.tsv | sed 's/\.tsv//')
do
    echo "==> ${group}"
    cat ${group}.tsv |
        tsv-filter --ne 5:0 |
        wc -l
    
    tsv-join all.mut.vcf -k 1,2,3 \
        --filter-file <(cat ${group}.tsv | tsv-filter --ne 5:0) |
        wc -l
done
```

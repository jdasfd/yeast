# Check wild yeast vcf

The [article](https://www.nature.com/articles/s41586-022-04823-w) contains thousands of mutations from 21 genes (synonymous or non-synonymous mutations). The goal is to validate all mutations among wild yeasts.

`mkdir -p ~/data/yeast/info`

Download all supplementary files and tables from the [website](https://www.nature.com/articles/s41586-022-04823-w) and save them into info dir.

Steps followed here are processes to convert the mutation info into vcf format for better analyze.

There are some data from the [README.md](https://github.com/wang-q/pars#readme). Better completing those steps.

Scripts could be partly found in the github repo [pars](https://github.com/wang-q/pars). Make sure you clone it before.

## Software

```bash
cd ~/Scripts
git clone https://github.com/wang-q/fig_table.git

cd fig_table
cp xlsx2csv.pl ~/data/yeast/scripts
```

```bash
brew install blast
brew install bcftools
brew install wang-q/tap/tsv-utils
brew install brewsci/bio/vt
```

## Random mutations in the article

### Extract 21 gene names and sequences

According to SGD, there is system name and standard name for one same gene. The standard name was adopted in the article. When it comes to different processes, two names were respectively useful in different aspects.

After that, the sequences were extracted from the mutation info.

- Gene names

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

# average fitness of muts were not all available
# remove those #DIV/0! in excel
perl scripts/xlsx2csv.pl -f info/41586_2022_4823_MOESM9_ESM.xlsx \
    --sheet "Fig. 2abc" |
    sed '1d' |
    tsv-select -d, -f 1,3,2 |
    mlr --icsv --otsv cat |
    tr '-' '\t' |
    tsv-filter --not-iregex 5:# \
    > info/fit.tsv

# count all mutations with available fitness data
cat info/fit.tsv | wc -l
#8341
```

- Gene sequences

```bash
# get fa according to the article info
# all detected muts were included
cd ~/data/yeast/gene

rm gene.mut.fa

for gene in $(cat stdname.lst)
do
    echo "==> ${gene}"

    echo ">${gene}" \
        >> gene.mut.fa

    cat ../info/fit.tsv |
        grep "^$gene" |
        tsv-select -f 1,2,3 |
        tsv-uniq |
        sort -nk 2,2 |
        perl -nae '
            if (!defined $i){
                $i = 1;
                print $F[2];
                $i++;
                }
            elsif($i == $F[1]){
                print $F[2];
                $i++;
            }
            else{
                print "N";
                $i++;
                redo;
            }
            END{print qq{\n};}' \
        >> gene.mut.fa
done

# All genes and their length
faops size gene.mut.fa
# check whether all undetected muts converted to N
```

### Use blast to get genome location of genes

```bash
# S288c database was used in pars
# make sure you have already completed the most part of pars
# otherwise using the followed command

# formatdb
# cd ~/data/mrna-structure/blast
# makeblastdb -dbtype nucl -in S288c.fa -parse_seqids

cd ~/data/yeast/gene

# blast every transcripts against genome
blastn -task blastn -evalue 1e-3 -num_threads 4 -num_descriptions 10 -num_alignments 10 -outfmt 0 \
    -dust yes -soft_masking true \
    -db ~/data/mrna-structure/blast/S288c.fa -query gene.mut.fa -out gene.blast

# identity 90 could give out the most gene wanted
perl ~/Scripts/pars/blastn_transcript.pl -i 70 -f gene.blast -m 0

cat gene.blast.tsv | wc -l
#19
# LSM1 and VMA7 excluded
# the reason was that 2 genes were not aligned in full length

# add those part manually according to gene.blast
echo -e "LSM1\t150\tX\t187438\t187587\t-" >> gene.blast.tsv
echo -e "VMA21\t150\tVII\t698635\t698784\t+" >> gene.blast.tsv
```

### Get random mutations of 21 genes in vcf format

There are few things should be considered:

1. The file format of 1002 genome project VCF is [VCF v4.1](http://samtools.github.io/hts-specs/VCFv4.1.pdf). In this format, all  info of mutations was recorded on the forward strand of reference yeast genome R64.

2. In the article, each position of mutations was recorded according to the coding region. So after I researched on the [genome browser](https://browse.yeastgenome.org/), it should be noticed that the blast result should be relocated on the forward strand.

In other words, for the forward strand (+), the genome location is alright, but it should be upside down for the reverse strand (-).

- Forward strand: left 5' -> right 3'.
- Reverse strand: right 5' -> left 3'.

```bash
cd ~/data/yeast

# a script converting results to vcf format
perl scripts/loc2vcf.pl -b gene/gene.blast.tsv \
    -t info/fit.tsv |
    tsv-sort -nk 2,2 \
    > vcf/random.snp.tsv

cat vcf/random.snp.tsv | wc -l
#8341
# the same number of info/fit.tsv

cat vcf/random.snp.tsv |
    s
```

## SNP in wild population

### All SNPs in the whole genome

All information were included in a tsv: chr, pos, ALT, REF, freq, num_ALT, num_all.

```bash
cd ~/data/yeast

# echo 1573791+82429+2147 | bc
# all equations could be validated by the command

# get all SNPs
# bcftools norm to split mutations within a position
# vt decompose_blocksub will automatically change them to 1 base
bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 8 |
    bcftools norm -m -both |
    vt decompose_blocksub - |
    bcftools view -V indels |
    bcftools query -f \
    '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' \
    -o vcf/all.snp.tsv

# multimuts snp
cat vcf/all.snp.tsv | tsv-summarize -g 1,2 --count | tsv-summarize -g 3 --count
#1       1573791
#2       82429
#3       2147

# genome absolute positions with mutation(s)
for i in {1..3}
do
    cat vcf/all.snp.tsv | tsv-uniq -f 1,2 -a $i | wc -l
done
#1658367 (1573791+82429+2147)
#84576 (82429+2147)
#2147

# the number of all SNPs
cat vcf/all.snp.tsv | wc -l
#1745090 (1658367+84576+2147 = 1573791+82429*2+2147*3)

# freq >= 0.05
cat vcf/all.snp.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       137704
#2       731
#3       5

for i in {1..3}
do
    cat vcf/all.snp.tsv |
        tsv-filter --ge 5:0.05 |
        tsv-uniq -f 1,2 -a $i |
        wc -l
done
#138440 (137704+731+5)
#736 (731+5)
#5

cat vcf/all.snp.tsv | tsv-filter --ge 5:0.05 | wc -l
#139181 (138440+736+5 = 137704+731*2+5*3)
```

### SNPs in gene regions

```bash
cd ~/data/yeast

# all gene regions added to 1-based bed format
cat ../mrna-structure/gene-filter/gene_list.csv |
    perl -nla -F, -e '
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
    $F[1] =~ /(.+)\([+-]\):(\d+)-(\d+)/;
    $chr = "chromosome" . $roman{$1};
    $begin = $2;
    $end = $3;
    print "$chr\t$begin\t$end";
    ' \
    > gene/gene_all.bed

cat gene/gene_all.bed | wc -l
#6579
# 6579 genes in total

bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 6 -R gene/gene_all.bed |
    bcftools norm -m -both |
    vt decompose_blocksub - |
    bcftools view -V indels |
    bcftools query -f \
    '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' \
    -o vcf/gene_all.snp.tsv

cat vcf/gene_all.snp.tsv |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       1027584
#2       45004
#3       917

for i in {1..3}
do
    cat vcf/gene_all.snp.tsv | tsv-uniq -f 1,2 -a $i | wc -l
done
#1073505 (1027584+45004+917)
#45921 (45004+917)
#917

cat vcf/gene_all.snp.tsv | wc -l
#1120343 (1073505+45921+917 = 1027584+45004*2+917*3)

# freq >= 0.05
cat vcf/gene_all.snp.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       82319
#2       311

for i in {1..3}
do
    cat vcf/gene_all.snp.tsv |
        tsv-filter --ge 5:0.05 |
        tsv-uniq -f 1,2 -a $i |
        wc -l
done
#82630 (82319+311+0)
#311 (311+0)
#0

cat vcf/gene_all.snp.tsv | tsv-filter --ge 5:0.05 | wc -l
#82941 (82630+311+0 = 82319+311*2+0*3)
```

### SNPs not in gene regions

```bash
cd ~/data/yeast

cat vcf/all.snp.tsv | tsv-join -k 1,2,3,4 -e \
    -f vcf/gene_all.snp.tsv \
    > vcf/not_gene.snp.tsv

cat vcf/not_gene.snp.tsv |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       546207
#2       37425
#3       1230

for i in {1..3}
do
    cat vcf/not_gene.snp.tsv | tsv-uniq -f 1,2 -a $i | wc -l
done
#584862 (546207+37425+1230)
#38655 (37425+1230)
#1230

cat vcf/not_gene.snp.tsv | wc -l
#624747 (1073505+45921+917 = 1027584+45004*2+917*3)

echo 1120343+624747 | bc
#1745090
# all.snp.tsv was seperated into 2 parts

# freq >= 0.05
cat vcf/not_gene.snp.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       55385
#2       420
#3       5

for i in {1..3}
do
    cat vcf/not_gene.snp.tsv |
        tsv-filter --ge 5:0.05 |
        tsv-uniq -f 1,2 -a $i |
        wc -l
done
#55810 (55385+420+5)
#425 (420+5)
#5

cat vcf/not_gene.snp.tsv | tsv-filter --ge 5:0.05 | wc -l
#56240 (55810+425+5 = 55385+420*2+5*3)
```

### SNPs in selected 21 genes

```bash
cd ~/data/yeast

# bcftools required 1-based bed for snp filtering
cat gene/gene.blast.tsv |
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
          my $chr = $roman{$F[2]};
          print "chromosome$chr\t$F[3]\t$F[4]";
     ' \
     > gene/gene_21.bed

bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 6 -R gene/gene_21.bed |
    bcftools norm -m -both |
    vt decompose_blocksub - |
    bcftools view -V indels |
    bcftools query -f \
    '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' \
    -o vcf/gene_21.snp.tsv

cat vcf/gene_21.snp.tsv |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       286
#2       10

for i in {1..3}
do
    cat vcf/gene_21.snp.tsv | tsv-uniq -f 1,2 -a $i | wc -l
done
#296 (286+10+0)
#10 (10+0)
#0

cat vcf/gene_21.snp.tsv | wc -l
#306 (296+10+0 = 286+10*2+0*3)

cat vcf/gene_21.snp.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       22

for i in {1..3}
do
    cat vcf/gene_21.snp.tsv |
        tsv-filter --ge 5:0.05 |
        tsv-uniq -f 1,2 -a $i |
        wc -l
done
#22
#0
#0

cat vcf/gene_21.snp.tsv | tsv-filter --ge 5:0.05 | wc -l
#22
```

### Chi-square tests

Check if there is any differences between regions.

```bash
bash ~/data/yeast/scripts/chi.sh
```

Result:

```txt
==> genome vs gene
         [,1]    [,2]
[1,] 12071326 8813103
[2,]  1658367 1073505

        Pearson's Chi-squared test with Yates' continuity correction

data:  x
X-squared = 8369.9, df = 1, p-value < 2.2e-16

==> gene vs other
        [,1]    [,2]
[1,] 8813103 3258223
[2,] 1073505  584862

        Pearson's Chi-squared test with Yates' continuity correction

data:  x
X-squared = 49545, df = 1, p-value < 2.2e-16

==> gene vs 21
        [,1] [,2]
[1,] 8813103 3138
[2,] 1073505  296

        Pearson's Chi-squared test with Yates' continuity correction

data:  x
X-squared = 17.542, df = 1, p-value = 2.81e-05

```

### Random mutations among wild groups

```bash
cd ~/data/yeast/vcf

cat random.snp.tsv |
    tsv-join -k 1,2,3,4 \
    -f gene_21.snp.tsv -a 5,6,7 \
    > random.wild.snp.tsv

cat random.snp.tsv |
    tsv-join -k 1,2,3,4 \
    -f gene_21.snp.tsv -e \
    > random.not_wild.snp.tsv

cat random.wild.snp.tsv |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       262
#2       10

cat random.wild.snp.tsv |
    tsv-filter --ge 8:0.05 |
    tsv-summarize -g 1,2 --count |
    tsv-summarize -g 3 --count
#1       20
```

### Statistical results of mutations

```bash
cd ~/data/yeast/vcf

cat random.wild.snp.tsv |
    tsv-summarize -g 6 --mean 5 --median 5 |
    sed '1itype\tfit_mean\tfit_median' |
    mlr --itsv --omd cat

cat random.wild.snp.tsv |
    tsv-select -f 6,7,5,8 |
    perl -nla -e '
        print join("\t",@F) if $F[0] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[0] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1itype\tgene\tfit\tfreq' \
    > ../results/wild.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    library(plyr)
    args <- commandArgs(T)
    wild <- read_tsv(args[1], show_col_types = FALSE)
    wildv <- ddply(wild, "type", summarise, grp.mean = mean(fit))
    p <- ggplot(wild, aes(x = fit, fill = type)) +
         geom_histogram(binwidth = 0.0025, alpha = 0.5, position = "identity") +
         geom_vline(data = wildv, aes(xintercept = grp.mean, color = type), linetype = "dashed")
    ggsave(p, height = 6, width = 6, file = "../results/wild.all.fit.pdf")
' ../results/wild.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    library(gridExtra)
    args <- commandArgs(T)
    wild <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(wild, aes(x = fit, fill = type)) +
         geom_histogram(alpha = 0.5, position = "identity") +
         facet_wrap(~gene, nrow = 6)
    ggsave(p, height = 15, width = 15, file = "../results/wild.gene.fit.pdf")
' ../results/wild.tsv
```

### 

## Split strains from 1002 genomes project into subpopulations

### Get wild yeasts source from 1002 genomes project

- Ecological origins

The groups were divided directly according to the info recorded in the project, which contained each strain's ecological origin. See the [project main page](http://1002genomes.u-strasbg.fr/) for more details.

```bash
mkdir ~/data/yeast/isolates
cd ~/data/yeast/isolates

wget http://1002genomes.u-strasbg.fr/isolates/page8/files/1002genomes.txt

# Remove references part
# Change the col2 to variables
cat 1002genomes.txt |
    sed '1013,1042d' |
    tsv-select -f 1,4 |
    sed 's/Human, clinical/Clinical/' |
    tr ' ' '_' \
    > 1002genomes.tsv

cat 1002genomes.tsv |
    sed '1d' |
    cut -f 2 |
    sort |
    uniq \
    > group.lst

for catgry in $(cat group.lst)
do
    echo "==> ${catgry}"
    cat 1002genomes.tsv |
        sed '1d' |
        tsv-filter --str-eq 2:${catgry} |
        tsv-select -f 1 \
        > ${catgry}.lst
done

cat 1002genomes.tsv | sed '1d' | wc -l
#1011

wc -l *.lst
#  23 group.lst
#1034 total
# the result is 1011 after subtraction, split correctly
```

- Lineages (update soon)

### Split the `1011Matrix.gvcf.gz` into small files according to subpopulations

```bash
cd ~/data/yeast
mkdir -p vcf/group

# split vcf according to groups
for group in $(cat isolates/group.lst)
do
    echo "==> ${group}"
    
    sample=$(cat isolates/${group}.lst |
        tr '\n' ',' |
        sed 's/,$//')
    
    bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz \
        --threads 8 -s ${sample} |
        bcftools +fill-tags -Ob -o vcf/group/1011Matrix.${group}.bcf
done
# -Ob: output bcf (compressed)
# only compressed bcf format could be indexed

cd ~/data/yeast/vcf/group

parallel -j 4 " \
    bcftools index --threads 3 {} \
" ::: $(ls *.bcf)
## Genome alignment

### Download genomes

```bash
mkdir -p ~/data/yeast/download
cd ~/data/yeast/download

# BY4742 strain was used in the article
wget http://sgd-archive.yeastgenome.org/sequence/strains/BY4742/BY4742_Toronto_2012/BY4742_Toronto_2012.fsa.gz
wget http://sgd-archive.yeastgenome.org/sequence/strains/BY4742/BY4742_Toronto_2012/BY4742_Toronto_2012.gff.gz

# S288c - reference genome
aria2c -c ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
aria2c -c ftp://ftp.ensembl.org/pub/release-105/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz
```

### Prepare sequences

```bash
cd ~/data/yeast

# reference
egaz prepseq \
    download/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz \
    --repeatmasker "--species Fungi --gff --parallel 12" \
    --min 1000 --gi -v \
    -o GENOMES/S288c

egaz prepseq \
    download/BY4742_Toronto_2012.fsa.gz \
    --repeatmasker "--species Fungi --gff --parallel 12" \
    --min 1000 --gi -v \
    -o GENOMES/BY4742

gzip -dcf download/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz > GENOMES/S288c/chr.gff
gzip -dcf download/BY4742_Toronto_2012.gff.gz > GENOMES/BY4742/chr.gff

# prep assembly
egaz template \
    download \
    --prep -o GENOMES \
    --min 1000 --about 1_000_000 \
    -v --repeatmasker "--species Fungi --parallel 12"

bash GENOMES/0_prep.sh
```

## Statistical analysis

### Basic info

- The numbers of SNPs occurred among subpopulations

```bash
mkdir ~/data/yeast/results
cd ~/data/yeast/vcf

cat all.vcf.tsv | wc -l
#819
# totally 819 snps

cat all.vcf.tsv | tsv-filter --ge 5:0.05 | wc -l
#431
# 431 snps population freq >= 0.05

# uniq all snps
cat all.vcf.tsv |
    tsv-select -f 1,2,3,4 |
    tsv-uniq |
    wc -l
#248 (8004-7756, right)
# totally 248 snps found among wild groups

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-select -f 1,2,3,4 |
    tsv-uniq |
    wc -l
#112
# totally 112 high freq (>= 0.05) snps found among groups

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 11 --count |
    tsv-join -f <(cat all.vcf.tsv | tsv-summarize -g 11 --count) -k 1 -a 2 |
    tsv-sort -nk 3,3 -r |
    sed '1igroup\thigh_freq\tall' |
    mlr --itsv --omd cat

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 11 --count |
    awk '{print ($0 "\thigh_freq")}' |
    sed '1igroup\tnum\tcatgry' > tmp &&
cat all.vcf.tsv |
    tsv-summarize -g 11 --count |
    awk '{print ($0 "\tall")}' >> tmp &&
mv tmp ../results/group.num.tsv

# plot
# need to change x axis
Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    data <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(data, aes(x = reorder(group, -num), y = num, fill = catgry)) +
         geom_bar(stat="identity", position=position_dodge()) +
         geom_text(aes(label = num), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
         scale_fill_brewer(palette="Paired") +
         theme(axis.text.x = element_text(angle = 315))
    ggsave(p, height = 6, width = 15, file = "../results/group.num.pdf")
' ../results/group.num.tsv
```

| group        | high_freq | all |
|--------------|-----------|-----|
| Nature       | 23        | 83  |
| Wine         | 10        | 63  |
| Beer         | 24        | 59  |
| Tree         | 22        | 57  |
| Fruit        | 21        | 54  |
| Clinical     | 20        | 54  |
| Unknown      | 17        | 39  |
| Distillery   | 18        | 38  |
| Palm_wine    | 22        | 37  |
| Insect       | 33        | 37  |
| Fermentation | 20        | 34  |
| Water        | 27        | 33  |
| Soil         | 27        | 33  |
| Industrial   | 23        | 32  |
| Human        | 16        | 29  |
| Bakery       | 22        | 29  |
| Bioethanol   | 12        | 20  |
| Flower       | 18        | 19  |
| Dairy        | 13        | 19  |
| Cider        | 17        | 19  |
| Sake         | 11        | 16  |
| Lab_strain   | 9         | 9   |
| Probiotic    | 6         | 6   |

- The numbers of subpopulations of an SNP

```bash
cd ~/data/yeast/vcf

cat ../isolates/group.lst | wc -l
#23

# all snps occurred at least 10 groups
cat all.vcf.tsv |
    tsv-summarize -g 1,2,3,4,10 --count |
    tsv-sort -r -nk 6,6 |
    tsv-select -f 1,5,6 |
    tsv-filter --ge 3:10 |
    sed '1ichr\tgene\tgroup' |
    mlr --itsv --omd cat
```

| chr          | gene  | group |
|--------------|-------|-------|
| chromosome12 | EST1  | 23    |
| chromosome7  | GET1  | 22    |
| chromosome7  | GET1  | 22    |
| chromosome12 | TSR2  | 22    |
| chromosome6  | RPL29 | 21    |
| chromosome5  | IES6  | 21    |
| chromosome14 | EOS1  | 21    |
| chromosome14 | EOS1  | 21    |
| chromosome3  | BUD23 | 20    |
| chromosome5  | IES6  | 19    |
| chromosome12 | CCW12 | 19    |
| chromosome12 | EST1  | 18    |
| chromosome3  | BUD23 | 16    |
| chromosome7  | GET1  | 13    |
| chromosome6  | RPL29 | 12    |
| chromosome5  | IES6  | 11    |
| chromosome13 | ASC1  | 11    |
| chromosome12 | EST1  | 11    |
| chromosome5  | IES6  | 10    |
| chromosome12 | TSR2  | 10    |

- Frequencies of each SNP occurred more than 10 subpopulations

```bash
mkdir ~/data/yeast/vcf/freq
cd ~/data/yeast/vcf

# freq among groups
cat all.vcf.tsv |
    tsv-select -f 1,2,3,4 |
    tsv-uniq |
    parallel --colsep '\t' -j 4 -k '
        cat all.vcf.tsv |
        tsv-filter --str-eq 1:{1} --eq 2:{2} --str-eq 3:{3} --str-eq 4:{4} |
        tsv-select -f 11,9,5 |
        sed "1igroup\ttype\tfreq" \
        > freq/{1}_{2}_{3}_{4}.tsv
    '

cd ~/data/yeast/vcf/freq

ls *.tsv | wc -l
#248

cat ../all.vcf.tsv |
    tsv-summarize -g 1,2,3,4 --count | 
    tsv-filter --ge 5:10 |
    wc -l
#20

for name in $(wc -l *.tsv |
              grep -v 'total$' |
              datamash reverse -W |
              tsv-filter --lt 2:11 | #header should be included
              tsv-select -f 1)
do
    rm ${name}
done

wc -l *.tsv |
    grep -v 'total$' |
    wc -l
#20

# plot
for file in $(ls *.tsv)
do
    echo "==>${file}"
    Rscript -e '
        library(ggplot2)
        library(readr)
        args <- commandArgs(T)
        freq <- read_tsv(args[1], show_col_types = FALSE)
        save <- paste0(args[1], ".pdf")
        p <- ggplot(freq, aes(x = reorder(group, -freq), y = freq)) +
             geom_bar(stat="identity") +
             geom_text(aes(label = freq), vjust=1.6, color="white",
                       position = position_dodge(0.9), size=3.5)+
             theme(axis.text.x = element_text(angle = 315))
        ggsave(p, height = 6, width = 10, file = save)
    ' ${file}
done
```

### Frequencies and fitness of SNP

- The mean and median of all mutation frequencies

```bash
cd ~/data/yeast/vcf

# freq mean and median
# snps existent in wild groups
cat all.vcf.tsv |
    tsv-summarize -g 9 --count --mean 5 --median 5 |
    sed '1imut_type\tnum\tfreq_mean\tfreq_median' |
    mlr --itsv --omd cat

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 9 --count --mean 5 --median 5 |
    sed '1imut_type\tnum\thigh_freq_mean\thigh_freq_median' |
    mlr --itsv --omd cat
```

| mut_type               | num | freq_mean      | freq_median |
|------------------------|-----|----------------|-------------|
| Synonymous_mutation    | 504 | 0.248842422897 | 0.0669872   |
| Nonsynonymous_mutation | 315 | 0.110011291778 | 0.0357143   |

| mut_type               | num | high_freq_mean | high_freq_median |
|------------------------|-----|----------------|------------------|
| Synonymous_mutation    | 292 | 0.41346027637  | 0.3074075        |
| Nonsynonymous_mutation | 139 | 0.221710007194 | 0.109375         |

```bash
cd ~/data/yeast/vcf

# output into a tsv for plot
cat all.vcf.tsv |
    tsv-summarize -g 11,9 --count --mean 5 --median 5 |
    sed '1igroup\ttype\tnum\tfreq_mean\tfreq_median' |
    mlr --itsv --omd cat

cat all.vcf.tsv |
    tsv-select -f 11,5,9 |
    perl -nla -e '
        print join("\t",@F) if $F[2] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[2] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1igroup\tfreq\ttype' \
    > ../results/group.freq.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    freq <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(freq, aes(x = type, y = freq, fill = type))+
         geom_boxplot() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    ggsave(p, height = 6, width = 15, file = "../results/group.freq.pdf")
' ../results/group.freq.tsv

cat all.vcf.tsv |
    tsv-summarize -g 11,9 --count --mean 5 --median 5 |
    perl -nla -e '
        print join("\t",@F) if $F[1] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[1] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1igroup\ttype\tnum\tfreq_mean\tfreq_median' \
    > ../results/group.freq.bar.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    freq <- read_tsv(args[1], show_col_types = FALSE)
    p1 <- ggplot(freq, aes(x = type, y = freq_mean, fill = type))+
         geom_col() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    p2 <- ggplot(freq, aes(x = type, y = freq_median, fill = type))+
         geom_col() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    ggsave(p1, height = 6, width = 15, file = "../results/group.freq.mean.bar.pdf")
    ggsave(p2, height = 6, width = 15, file = "../results/group.freq.median.bar.pdf")
' ../results/group.freq.bar.tsv

# >= 0.05
cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 11,9 --count --mean 5 --median 5 |
    sed '1igroup\ttype\tnum\thighfreq_mean\thighfreq_median' |
    mlr --itsv --omd cat

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-select -f 11,5,9 |
    perl -nla -e '
        print join("\t",@F) if $F[2] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[2] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1igroup\tfreq\ttype' \
    > ../results/group.highfreq.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    freq <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(freq, aes(x = type, y = freq, fill = type))+
         geom_boxplot() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    ggsave(p, height = 6, width = 15, file = "../results/group.highfreq.pdf")
' ../results/group.highfreq.tsv

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 11,9 --count --mean 5 --median 5 |
    perl -nla -e '
        print join("\t",@F) if $F[1] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[1] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1igroup\ttype\tnum\thighfreq_mean\thighfreq_median' \
    > ../results/group.highfreq.bar.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    freq <- read_tsv(args[1], show_col_types = FALSE)
    p1 <- ggplot(freq, aes(x = type, y = highfreq_mean, fill = type))+
         geom_col() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    p2 <- ggplot(freq, aes(x = type, y = highfreq_median, fill = type))+
         geom_col() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    ggsave(p1, height = 6, width = 15, file = "../results/group.highfreq.mean.bar.pdf")
    ggsave(p2, height = 6, width = 15, file = "../results/group.highfreq.median.bar.pdf")
' ../results/group.highfreq.bar.tsv

# distribution of freq
# all unique detected snps
cat all.vcf.tsv |
    tsv-select -f 11,5,9 |
    sed '1igroup\tfreq\ttype' \
    > ../results/dist.freq.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    library(plyr)
    library(gridExtra)
    args <- commandArgs(T)
    freq <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(freq, aes(x = freq, fill = type)) +
         geom_histogram(alpha = 0.5, position = "identity") +
         facet_wrap(~group, nrow = 6)
    ggsave(p, height = 15, width = 15, file = "../results/dist.freq.pdf")
' ../results/dist.freq.tsv
```

| group        | type                   | num | freq_mean       | freq_median |
|--------------|------------------------|-----|-----------------|-------------|
| Bakery       | Synonymous_mutation    | 15  | 0.31441444      | 0.27027     |
| Bakery       | Nonsynonymous_mutation | 14  | 0.102757864286  | 0.0675676   |
| Beer         | Synonymous_mutation    | 33  | 0.175525345455  | 0.0423729   |
| Beer         | Nonsynonymous_mutation | 26  | 0.0634637838462 | 0.0338983   |
| Bioethanol   | Synonymous_mutation    | 12  | 0.430555341667  | 0.304131    |
| Bioethanol   | Nonsynonymous_mutation | 8   | 0.136574075     | 0.03703705  |
| Cider        | Nonsynonymous_mutation | 9   | 0.216298922222  | 0.205882    |
| Cider        | Synonymous_mutation    | 10  | 0.38823524      | 0.2058825   |
| Clinical     | Nonsynonymous_mutation | 23  | 0.0739591773913 | 0.0233645   |
| Clinical     | Synonymous_mutation    | 31  | 0.146980374839  | 0.0280374   |
| Dairy        | Nonsynonymous_mutation | 6   | 0.370370266667  | 0.1944443   |
| Dairy        | Synonymous_mutation    | 13  | 0.410826238462  | 0.111111    |
| Distillery   | Nonsynonymous_mutation | 13  | 0.119363384615  | 0.0344828   |
| Distillery   | Synonymous_mutation    | 25  | 0.2461576       | 0.0517241   |
| Fermentation | Synonymous_mutation    | 22  | 0.332702027273  | 0.277778    |
| Fermentation | Nonsynonymous_mutation | 12  | 0.15625         | 0.0277778   |
| Flower       | Synonymous_mutation    | 13  | 0.368131892308  | 0.214286    |
| Flower       | Nonsynonymous_mutation | 6   | 0.202380916667  | 0.0892858   |
| Fruit        | Synonymous_mutation    | 36  | 0.150369286111  | 0.0425532   |
| Fruit        | Nonsynonymous_mutation | 18  | 0.0721040277778 | 0.0319149   |
| Human        | Nonsynonymous_mutation | 14  | 0.221198171429  | 0.02419355  |
| Human        | Synonymous_mutation    | 15  | 0.43870974      | 0.193548    |
| Industrial   | Nonsynonymous_mutation | 14  | 0.0955782714286 | 0.06904765  |
| Industrial   | Synonymous_mutation    | 18  | 0.274338588889  | 0.2         |
| Insect       | Nonsynonymous_mutation | 9   | 0.116666666667  | 0.05        |
| Insect       | Synonymous_mutation    | 28  | 0.227678571429  | 0.1         |
| Lab_strain   | Synonymous_mutation    | 8   | 0.625           | 0.5         |
| Lab_strain   | Nonsynonymous_mutation | 1   | 0.5             | 0.5         |
| Nature       | Synonymous_mutation    | 53  | 0.0990637173585 | 0.0192308   |
| Nature       | Nonsynonymous_mutation | 30  | 0.0497360686667 | 0.0192308   |
| Palm_wine    | Synonymous_mutation    | 26  | 0.290384630769  | 0.0666667   |
| Palm_wine    | Nonsynonymous_mutation | 11  | 0.148484872727  | 0.0333333   |
| Probiotic    | Nonsynonymous_mutation | 3   | 0.5             | 0.5         |
| Probiotic    | Synonymous_mutation    | 3   | 1               | 1           |
| Sake         | Synonymous_mutation    | 13  | 0.688216115385  | 0.989362    |
| Sake         | Nonsynonymous_mutation | 3   | 0.386524866667  | 0.148936    |
| Soil         | Synonymous_mutation    | 20  | 0.255921065     | 0.13815815  |
| Soil         | Nonsynonymous_mutation | 13  | 0.106275330769  | 0.0789474   |
| Tree         | Synonymous_mutation    | 41  | 0.153325256098  | 0.03125     |
| Tree         | Nonsynonymous_mutation | 16  | 0.08081440625   | 0.03125     |
| Unknown      | Nonsynonymous_mutation | 19  | 0.0808270684211 | 0.0357143   |
| Unknown      | Synonymous_mutation    | 20  | 0.251785705     | 0.0714286   |
| Water        | Nonsynonymous_mutation | 15  | 0.09298244      | 0.0526316   |
| Water        | Synonymous_mutation    | 18  | 0.269736877778  | 0.0789473   |
| Wine         | Synonymous_mutation    | 31  | 0.108155583871  | 0.00806452  |
| Wine         | Nonsynonymous_mutation | 32  | 0.0402662403125 | 0.007056455 |

| group        | type                   | num | highfreq_mean  | highfreq_median |
|--------------|------------------------|-----|----------------|-----------------|
| Bakery       | Synonymous_mutation    | 11  | 0.420147463636 | 0.364865        |
| Bakery       | Nonsynonymous_mutation | 11  | 0.127097236364 | 0.0810811       |
| Beer         | Nonsynonymous_mutation | 9   | 0.1447334      | 0.059322        |
| Beer         | Synonymous_mutation    | 15  | 0.357907153333 | 0.245763        |
| Bioethanol   | Nonsynonymous_mutation | 4   | 0.25462965     | 0.2592595       |
| Bioethanol   | Synonymous_mutation    | 8   | 0.6319441375   | 0.7129625       |
| Cider        | Nonsynonymous_mutation | 9   | 0.216298922222 | 0.205882        |
| Cider        | Synonymous_mutation    | 8   | 0.4779411      | 0.3088235       |
| Clinical     | Nonsynonymous_mutation | 7   | 0.202879742857 | 0.122642        |
| Clinical     | Synonymous_mutation    | 13  | 0.333597338462 | 0.264423        |
| Dairy        | Nonsynonymous_mutation | 4   | 0.5370369      | 0.583333        |
| Dairy        | Synonymous_mutation    | 9   | 0.576954788889 | 0.740741        |
| Distillery   | Nonsynonymous_mutation | 5   | 0.2586206      | 0.241379        |
| Distillery   | Synonymous_mutation    | 13  | 0.446854769231 | 0.517241        |
| Fermentation | Synonymous_mutation    | 15  | 0.474074073333 | 0.472222        |
| Fermentation | Nonsynonymous_mutation | 5   | 0.34166664     | 0.25            |
| Flower       | Synonymous_mutation    | 13  | 0.368131892308 | 0.214286        |
| Flower       | Nonsynonymous_mutation | 5   | 0.23571424     | 0.107143        |
| Fruit        | Nonsynonymous_mutation | 5   | 0.19148938     | 0.106383        |
| Fruit        | Synonymous_mutation    | 16  | 0.30563418125  | 0.209528        |
| Human        | Nonsynonymous_mutation | 6   | 0.49193555     | 0.564516        |
| Human        | Synonymous_mutation    | 10  | 0.6467743      | 0.806452        |
| Industrial   | Nonsynonymous_mutation | 9   | 0.130158811111 | 0.1             |
| Industrial   | Synonymous_mutation    | 14  | 0.344387714286 | 0.2416665       |
| Insect       | Nonsynonymous_mutation | 7   | 0.142857142857 | 0.05            |
| Insect       | Synonymous_mutation    | 26  | 0.243269230769 | 0.125           |
| Lab_strain   | Synonymous_mutation    | 8   | 0.625          | 0.5             |
| Lab_strain   | Nonsynonymous_mutation | 1   | 0.5            | 0.5             |
| Nature       | Nonsynonymous_mutation | 7   | 0.145846857143 | 0.0882353       |
| Nature       | Synonymous_mutation    | 16  | 0.28485573125  | 0.216346        |
| Palm_wine    | Nonsynonymous_mutation | 5   | 0.29000008     | 0.166667        |
| Palm_wine    | Synonymous_mutation    | 17  | 0.426470629412 | 0.183333        |
| Probiotic    | Nonsynonymous_mutation | 3   | 0.5            | 0.5             |
| Probiotic    | Synonymous_mutation    | 3   | 1              | 1               |
| Sake         | Synonymous_mutation    | 9   | 0.988179777778 | 1               |
| Sake         | Nonsynonymous_mutation | 2   | 0.569149       | 0.569149        |
| Soil         | Synonymous_mutation    | 19  | 0.268005552632 | 0.184211        |
| Soil         | Nonsynonymous_mutation | 8   | 0.157894775    | 0.0789474       |
| Tree         | Synonymous_mutation    | 18  | 0.318424888889 | 0.230469        |
| Tree         | Nonsynonymous_mutation | 4   | 0.250992       | 0.09375         |
| Unknown      | Nonsynonymous_mutation | 6   | 0.1904762      | 0.142857        |
| Unknown      | Synonymous_mutation    | 11  | 0.435064927273 | 0.428571        |
| Water        | Synonymous_mutation    | 16  | 0.3001645125   | 0.131579        |
| Water        | Nonsynonymous_mutation | 11  | 0.117224854545 | 0.0789474       |
| Wine         | Nonsynonymous_mutation | 6   | 0.169054333333 | 0.09173405      |
| Wine         | Synonymous_mutation    | 4   | 0.75201625     | 0.947581        |

- The mean and median of all mutation fitness

```bash
cd ~/data/yeast/vcf

# fitness mean and median
# uniq snps from different groups
# snps existent in wild groups
cat all.vcf.tsv |
    tsv-uniq -f 1,2,3,4 |
    tsv-summarize -g 9 --count --mean 8 --median 8 |
    sed '1itype\tnum\tfit_mean\tfit_median' |
    mlr --itsv --omd cat

# snps nonexistent in wild groups
cat other.vcf.tsv |
    tsv-summarize -g 6 --count --mean 5 --median 5 |
    sed '1itype\tnum\tfit_mean\tfit_median' |
    mlr --itsv --omd cat

# distribution of fitness
# all unique detected snps
cat all.vcf.tsv |
    tsv-uniq -f 1,2,3,4 |
    tsv-select -f 8,9 |
    sed '1ifit\ttype' \
    > ../results/dist.fit.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    library(plyr)
    args <- commandArgs(T)
    fit <- read_tsv(args[1], show_col_types = FALSE)
    fitv <- ddply(fit, "type", summarise, grp.mean = mean(fit))
    p <- ggplot(fit, aes(x = fit, fill = type)) +
         geom_histogram(alpha = 0.5, position = "identity") +
         geom_vline(data = fitv, aes(xintercept = grp.mean, color = type), linetype = "dashed")
    ggsave(p, height = 6, width = 12, file = "../results/dist.fit.pdf")
' ../results/dist.fit.tsv

# plot
cat all.vcf.tsv |
    tsv-uniq -f 1,2,3,4 |
    tsv-select -f 8,9 |
    awk '{print ("yes\t" $0)}' |
    perl -nla -e '
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Non.+$/N_mut/;
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Sy.+$/S_mut/;
        ' |
    sed '1iexist\tfit\ttype' \
    > ../results/all.fit.tsv
 
cat other.vcf.tsv |
    tsv-filter --str-ne 6:Nonsense_mutation |
    tsv-select -f 5,6 |
    awk '{print ("no\t" $0)}' |
    perl -nla -e '
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Non.+$/N_mut/;
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Sy.+$/S_mut/;
        ' \
    >> ../results/all.fit.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1], show_col_types = FALSE)
p <- ggplot(fit, aes(x = type, y = fit, fill = type)) +
     geom_boxplot() +
     facet_grid(~exist) +
     ylim(0.96, 1.02) +
     theme(axis.text.x = element_text(angle = 315))
ggsave(p, height = 6, width = 15, file = "../results/all.fit.pdf")
' ../results/all.fit.tsv
```

| type                   | num | fit_mean       | fit_median     |
|------------------------|-----|----------------|----------------|
| Synonymous_mutation    | 136 | 0.988034235352 | 0.98874490225  |
| Nonsynonymous_mutation | 112 | 0.988007605256 | 0.988724384625 |


| type                   | num  | fit_mean       | fit_median   |
|------------------------|------|----------------|--------------|
| Nonsynonymous_mutation | 5929 | 0.984781866716 | 0.9879315165 |
| Synonymous_mutation    | 1665 | 0.98783557355  | 0.9887744535 |
| Nonsense_mutation      | 162  | 0.935173157813 | 0.9403812345 |

```bash
cd ~/data/yeast/vcf

cat other.vcf.tsv |
    tsv-filter --str-ne 6:Nonsense_mutation |
    tsv-select -f 5,6 |
    awk '{print ("other\t" $0)}' \
    > ../results/group.fit.tsv

cat all.vcf.tsv |
    tsv-select -f 11,8,9 \
    >> ../results/group.fit.tsv

cat ../results/group.fit.tsv |
    perl -nla -e '
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Non.+$/N_mut/;
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1igroup\tfit\ttype' \
    > tmp && \
    mv tmp ../results/group.fit.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1], show_col_types = FALSE)
p <- ggplot(fit, aes(x = type, y = fit, fill = type)) +
     geom_boxplot() +
     facet_grid(~group) +
     ylim(0.96, 1.02) +
     theme(axis.text.x = element_text(angle = 315))
ggsave(p, height = 6, width = 15, file = "../results/group.fit.pdf")
' ../results/group.fit.tsv

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-select -f 11,8,9 |
    perl -nla -e '
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Non.+$/N_mut/;
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Sy.+$/S_mut/;
        ' |
    sed '1igroup\tfit\ttype' \
    > ../results/group.highfreq.fit.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1], show_col_types = FALSE)
p <- ggplot(fit, aes(x = type, y = fit, fill = type)) +
     geom_boxplot() +
     facet_grid(~group) +
     theme(axis.text.x = element_text(angle = 315))
ggsave(p, height = 6, width = 15, file = "../results/group.highfreq.fit.pdf")
' ../results/group.highfreq.fit.tsv
```

### Chi-square

The experiment from the original article was almost a simulation of random mutations. All detected mutations were not biased to either nonsynonymous mutation or synonymous mutation. So all unique fixed SNPs after would be random as well if there were no other reasons.

```bash
cd ~/data/yeast/vcf

raw1=$(cat all.vcf.tsv |
          tsv-uniq -f 1,2,3,4 |
          tsv-summarize -g 9 --count |
          tsv-sort -nk 2,2 |
          datamash transpose |
          sed 1d)

raw2=$(cat other.vcf.tsv |
          tsv-summarize -g 6 --count |
          tsv-filter --str-ne 1:Nonsense_mutation |
          datamash transpose |
          sed 1d)

echo -e "$raw1\t$raw2" |
    parallel --colsep '\t' -j 1 -k '
        Rscript -e "
            x <- matrix(c({1},{3},{2},{4}), ncol = 2)
            old.warn <- options()$warn
            options(warn = -1)
            x
            chisq.test(x)
        "
'
#     [,1] [,2]
#[1,]  112  136
#[2,] 5929 1665

#        Pearson's Chi-squared test with Yates' continuity correction

#data:  x
#X-squared = 145.2, df = 1, p-value < 2.2e-16
```

### Linear regression

```bash
cd ~/data/yeast/vcf

cat all.vcf.tsv |
    tsv-select -f 11,5,8 |
    sed '1igroup\tfreq\tfit' \
    > ../results/group.freq_fit.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
data <- read_tsv(args[1], show_col_types = FALSE)
p <- ggplot(data, aes(x = freq, y = fit)) +
     geom_point() +
     geom_smooth(method = "lm")
ggsave(p, height = 6, width = 8, file = "../results/group.freq_fit.mixed.pdf")
' ../results/group.freq_fit.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
data <- read_tsv(args[1], show_col_types = FALSE)
p <- ggplot(data, aes(x = freq, y = fit, color = group)) +
     geom_point() +
     geom_smooth(method = "lm") +
     facet_wrap(~group)
ggsave(p, height = 6, width = 15, file = "../results/group.freq_fit.pdf")
' ../results/group.freq_fit.tsv
```

## Count SNPs existed in other closely related species from original article

```bash
cd ~/data/yeast

# original data from the fig2e provided all muts found in other 5 yeast strains 
perl scripts/xlsx2csv.pl -f info/41586_2022_4823_MOESM9_ESM.xlsx \
    --sheet "Fig. 2e" |
    sed 1,3d |
    tsv-select -d, -f 3,5 |
    tsv-summarize -d, -g 1,2 --count |
    mlr --icsv --omd cat
```

| Nonsynonymous_mutation | No  | 5839 |
|------------------------|-----|------|
| Nonsynonymous_mutation | Yes | 169  |
| Synonymous_mutation    | No  | 1087 |
| Synonymous_mutation    | Yes | 714  |
| Nonsense_mutation      | No  | 146  |

"Yes" here meant muts could be found in at least one of nearest 5 yeast species. There were more than 800 muts existed among other yeast groups.

```bash
# test
bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 8 -V indels | bcftools norm -m- | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' | wc -l
Lines   total/split/realigned/skipped:  1625809/81320/0/0
#1709097

bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 8 -V indels | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' | wc -l
#1625809

bcftools view ../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 8 -V indels | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' | tsv-filter --iregex 4:, | wc -l
#81320

bcftools view ../../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 8 | bcftools norm -m- | bcftools view --threads 8 -V indels | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' | wc -l
Lines   total/split/realigned/skipped:  1754866/129507/0/0
#1745090

bcftools view ../../mrna-structure/vcf/1011Matrix.gvcf.gz -Ov --threads 8 | bcftools norm -m -both | vt decompose_blocksub - | bcftools view -V indels | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' | wc -l

#Lines   total/split/realigned/skipped:  1754866/129507/0/0
#
#stats: no. variants                       : 1920571
#       no. biallelic block substitutions  : 18744
#
#       no. additional SNPs                : 18744
#       no. variants after decomposition   : 1920571
#
#Time elapsed: 8m 11s
#
#1745090
```

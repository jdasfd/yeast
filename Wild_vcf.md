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
cat info/fit.tsv | tsv-uniq -f 1,2,3,4 > tmp && mv tmp info/fit.tsv
# check whether repeated

cat info/fit.tsv | wc -l
#8340
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
#8340
# the same number of info/fit.tsv

# diff fit
cat vcf/random.snp.tsv | tsv-filter --ge 5:1.01 > vcf/random.up.tsv
cat vcf/random.snp.tsv | tsv-filter --le 5:0.99 > vcf/random.down.tsv
cat vcf/random.snp.tsv | tsv-filter --gt 5:0.99 --lt 5:1.01 > vcf/random.neither.tsv

for change in {up,down,neither}
do
    echo "==> fitness ${change}"
    cat vcf/random.${change}.tsv |
        tsv-summarize -g 6 --count
done

#==> fitness up
#Synonymous_mutation     15
#Nonsynonymous_mutation  61
#==> fitness down
#Nonsynonymous_mutation  3669
#Synonymous_mutation     1043
#Nonsense_mutation       143
#==> fitness neither
#Nonsynonymous_mutation  2575
#Synonymous_mutation     808
#Nonsense_mutation       26
```

## SNP in wild population

### All SNPs in the whole genome

All information were included in a tsv: chr, pos, ALT, REF, freq, num_ALT, num_all.

```bash
cd ~/data/yeast

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

bash scripts/vcf_num.sh vcf/all.snp.tsv

# freq >= 0.5
cat vcf/all.snp.tsv | tsv-filter --ge 5:0.05 > vcf/all.high.snp.tsv

bash scripts/vcf_num.sh vcf/all.high.snp.tsv
```

| Mut        | Num     |
|------------|---------|
| All        | 1745090 |
| Pos        | 1658367 |
| One_SNP    | 1573791 |
| Two_SNPs   | 82429   |
| Three_SNPs | 2147    |

SNP that freq >= 0.05:

| Mut        | Num    |
|------------|--------|
| All        | 139181 |
| Pos        | 138440 |
| One_SNP    | 137704 |
| Two_SNPs   | 731    |
| Three_SNPs | 5      |

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

bash scripts/vcf_num.sh vcf/gene_all.snp.tsv

# freq >= 0.5
cat vcf/gene_all.snp.tsv | tsv-filter --ge 5:0.05 > vcf/gene_all.high.snp.tsv

bash scripts/vcf_num.sh vcf/gene_all.high.snp.tsv
```

| Mut        | Num     |
|------------|---------|
| All        | 1120343 |
| Pos        | 1073505 |
| One_SNP    | 1027584 |
| Two_SNPs   | 45004   |
| Three_SNPs | 917     |

SNP that freq >= 0.05:

| Mut        | Num   |
|------------|-------|
| All        | 82941 |
| Pos        | 82630 |
| One_SNP    | 82319 |
| Two_SNPs   | 311   |
| Three_SNPs |       |

### SNPs not in gene regions

```bash
cd ~/data/yeast

cat vcf/all.snp.tsv | tsv-join -k 1,2,3,4 -e \
    -f vcf/gene_all.snp.tsv \
    > vcf/not_gene.snp.tsv

bash scripts/vcf_num.sh vcf/not_gene.snp.tsv

# freq >= 0.5
cat vcf/not_gene.snp.tsv | tsv-filter --ge 5:0.05 > vcf/not_gene.high.snp.tsv

bash scripts/vcf_num.sh vcf/not_gene.high.snp.tsv
```

| Mut        | Num    |
|------------|--------|
| All        | 624747 |
| Pos        | 584862 |
| One_SNP    | 546207 |
| Two_SNPs   | 37425  |
| Three_SNPs | 1230   |

SNP that freq >= 0.05:

| Mut        | Num   |
|------------|-------|
| All        | 56240 |
| Pos        | 55810 |
| One_SNP    | 55385 |
| Two_SNPs   | 420   |
| Three_SNPs | 5     |

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

bash scripts/vcf_num.sh vcf/gene_21.snp.tsv

# freq >= 0.5
cat vcf/gene_21.snp.tsv | tsv-filter --ge 5:0.05 > vcf/gene_21.high.snp.tsv

bash scripts/vcf_num.sh vcf/gene_21.high.snp.tsv
```

| Mut        | Num |
|------------|-----|
| All        | 306 |
| Pos        | 296 |
| One_SNP    | 286 |
| Two_SNPs   | 10  |
| Three_SNPs |     |

SNP that freq >= 0.05:

| Mut        | Num |
|------------|-----|
| All        | 22  |
| Pos        | 22  |
| One_SNP    | 22  |
| Two_SNPs   |     |
| Three_SNPs |     |

### Chi-square tests

Check if there is any differences between regions.

```bash
bash ~/data/yeast/scripts/chi.sh
```

Result on the screen:

```txt
==> genome vs gene
        [,1]      [,2]
[1,] 1073505 0.7300857
[2,]  584862 0.2699143
[1] 0.6473266
[1] 0.3526734

        Chi-squared test for given probabilities

data:  x
X-squared = 57639, df = 1, p-value < 2.2e-16

==> gene vs 21
        [,1]         [,2]
[1,]     296 0.0003560607
[2,] 1073209 0.9996439393
[1] 0.0002757323
[1] 0.9997243

        Chi-squared test for given probabilities

data:  x
X-squared = 19.461, df = 1, p-value = 1.027e-05

```

### Random mutations among wild groups

```bash
cd ~/data/yeast/vcf

# muts exist or not in the wild group
cat random.snp.tsv |
    tsv-join -k 1,2,3,4 \
    -f gene_21.snp.tsv -a 5,6,7 \
    > random.wild.snp.tsv

cat random.snp.tsv |
    tsv-join -k 1,2,3,4 \
    -f gene_21.snp.tsv -e \
    > random.not_wild.snp.tsv

bash ../scripts/vcf_num.sh random.wild.snp.tsv

bash ../scripts/vcf_num.sh random.not_wild.snp.tsv
```

in the wild:

| Mut        | Num |
|------------|-----|
| All        | 282 |
| Pos        | 272 |
| One_SNP    | 262 |
| Two_SNPs   | 10  |
| Three_SNPs |     |

not in the wild:

| Mut        | Num  |
|------------|------|
| All        | 8058 |
| Pos        | 3047 |
| One_SNP    | 227  |
| Two_SNPs   | 629  |
| Three_SNPs | 2191 |

### Statistical results of mutations

There are some levels should be discussed:

1. Whether mutations exist among wild populations or not.
2. Whether mutations are synonymous or non-synonymous.

- Count and format table

```bash
cd ~/data/yeast/vcf

# muts in wild
cat random.wild.snp.tsv |
    tsv-summarize -g 6 --mean 5 --median 5 --count |
    sed '1itype\tfit_mean\tfit_median\tcount' |
    mlr --itsv --omd cat

# muts in wild (freq >= 0.05)
cat random.wild.snp.tsv |
    tsv-filter --ge 8:0.05 |
    tsv-summarize -g 6 --mean 5 --median 5 --count |
    sed '1itype\tfit_mean\tfit_median\tcount' |
    mlr --itsv --omd cat

# muts count (freq >= 0.1)
cat random.wild.snp.tsv |
    tsv-filter --ge 8:0.1 |
    tsv-summarize -g 6 --count |
    sed '1itype\tcount' |
    mlr --itsv --omd cat

cat random.not_wild.snp.tsv |
    tsv-summarize -g 6 --mean 5 --median 5 --count |
    sed '1itype\tfit_mean\tfit_median\tcount' |
    mlr --itsv --omd cat
```

Mutations fitness in wild:

| type                   | fit_mean       | fit_median    | count |
|------------------------|----------------|---------------|-------|
| Nonsynonymous_mutation | 0.988545219123 | 0.98948772825 | 130   |
| Synonymous_mutation    | 0.987720205444 | 0.988143842   | 152   |

High freq (>= 0.05):

| type                   | fit_mean       | fit_median    | count |
|------------------------|----------------|---------------|-------|
| Synonymous_mutation    | 0.985455825819 | 0.98603116925 | 13    |
| Nonsynonymous_mutation | 0.988684188175 | 0.990983907   | 7     |

Freq > 0.1 count:

| type                   | count |
|------------------------|-------|
| Synonymous_mutation    | 9     |
| Nonsynonymous_mutation | 2     |

Not in wild:

| type                   | fit_mean       | fit_median     | count |
|------------------------|----------------|----------------|-------|
| Nonsynonymous_mutation | 0.984878530688 | 0.988017172    | 6175  |
| Synonymous_mutation    | 0.987820844726 | 0.988834603375 | 1714  |
| Nonsense_mutation      | 0.933970180448 | 0.93983089275  | 169   |

- Plot them

```bash
mkdir ~/data/yeast/results
cd ~/data/yeast/vcf

# all random mutations exist among wild
cat random.wild.snp.tsv |
    tsv-select -f 6,7,5,8 |
    awk '{print $0 "\twild"}' |
    perl -nla -e '
        print join("\t",@F) if $F[0] =~ s/^Nonsy.+$/N_mut/;
        print join("\t",@F) if $F[0] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1itype\tgene\tfit\tfreq\texist' \
    > ../results/all.tsv

# other mutations
cat random.not_wild.snp.tsv |
    tsv-select -f 6,7,5 |
    awk '{print $0 "\t0\tnot"}' |
    perl -nla -e '
        next if $F[0] =~ /^Nonsense.+/;
        print join("\t",@F) if $F[0] =~ s/^Nonsy.+$/N_mut/;
        print join("\t",@F) if $F[0] =~ s/^Sy.+$/S_mut/;
    ' \
    >> ../results/all.tsv

# remove all nonsense mutations
cat ../results/all.tsv | wc -l
#8172
# echo 282+8058+1-169 | bc

# plot fit script
# mut type
Rscript ../scripts/dis_fit.r \
    -f ../results/all.tsv -t type -b 0.0025 \
    -o ../results/wild_type.fit.pdf

# exist or not
Rscript ../scripts/dis_fit.r \
    -f ../results/all.tsv -t exist -b 0.0025 \
    -o ../results/exist.fit.pdf

# 2 factors
Rscript ../scripts/dis_fit.r \
    -f ../results/all.tsv -t exist -b 0.0025 -s \
    -o ../results/exist.split.fit.pdf

# muts exist among wild with freq
cat random.wild.snp.tsv |
    tsv-select -f 6,7,5,8 |
    perl -nla -e '
        print join("\t",@F) if $F[0] =~ s/^Nonsy.+$/N_mut/;
        print join("\t",@F) if $F[0] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1itype\tgene\tfit\tfreq' \
    > ../results/wild.tsv

Rscript ../scripts/dis_freq.r \
    -f ../results/wild.tsv -b 0.05 \
    -o ../results/exist.freq.pdf

echo -e "type\tfit\tfreq" > ../results/change.tsv

for change in {up,down,neither}
do
    echo "==> fitness ${change}"
    cat random.wild.snp.tsv |
        tsv-join -k 1,2,3,4 -f random.${change}.tsv |
        tsv-select -f 5,8 |
        awk -v var=${change} '{print (var "\t" $0)}' \
        >> ../results/change.tsv
done

Rscript -e '
    library(ggplot2)
    library(readr)
    library(plyr)
    args <- commandArgs(T)
    wild <- read_tsv(args[1], show_col_types = FALSE)
    wildv <- ddply(wild, "type", summarise, grp.median = median(freq))
    p <- ggplot(wild, aes(x = freq, fill = type)) +
         geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity") +
         geom_vline(data = wildv, aes(xintercept = grp.median, color = type), linetype = "dashed")
    plm <- ggplot(wild, aes(x = freq, y = fit)) +
           geom_point() +
           geom_smooth(method = "lm")
    ggsave(p, height = 6, width = 6, file = "../results/change.freq.pdf")
    ggsave(plm, height = 6, width = 6, file = "../results/change.freq.lm.pdf")
' ../results/change.tsv
```

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

### Split the `1011Matrix.gvcf.gz` into small files according to subpopulations

```bash
cd ~/data/yeast
mkdir -p vcf/group

# number of strains
echo -e "group\tnum" > isolates/strain.tsv

for group in $(cat isolates/group.lst)
do
    echo "==> ${group} strain numbers"
    cat isolates/${group}.lst |
        wc -l |
        awk -v group=${group} '{print (group "\t" $0)}' \
        >> isolates/strain.tsv
done

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

for group in $(cat ../../isolates/group.lst)
do
    echo "==> ${group}"

    bcftools view 1011Matrix.${group}.bcf -Ov \
        --threads 6 -R ../../gene/gene_21.bed |
        bcftools norm -m -both |
        vt decompose_blocksub - |
        bcftools view -V indels |
        bcftools query -f \
        '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' | 
        tsv-filter --ne 7:0 |
        tsv-filter --ne 5:0 \
        > ../${group}.snp.tsv
done
```

### Statistical analysis

- Subpopulations basic info

```bash
cd ~/data/yeast/vcf

echo -e "Mut\tAll\tPos\tOne_SNP\tTwo_SNPs\tThree_SNPs" |
    datamash transpose > group.num.tsv

# combine all muts into one matrix-like tsv
for group in $(cat ../isolates/group.lst)
do
    echo "==> ${group}"
    cat group.num.tsv |
        tsv-join -H -k Mut \
        -f <(bash ../scripts/vcf_num.sh ${group}.snp.tsv |
                grep '^|' |
                grep -v '-' |
                sed 's/|  |/| 0 |/' |
                tsv-select -d '|' -f 2,3 |
                sed 's/^\s//' |
                sed 's/\s$//' |
                sed 's/ | /\t/' |
                sed "s/Num/${group}/") \
        -a ${group} > tmp && mv tmp group.num.tsv
done

# count basic mutations among groups
cat group.num.tsv | mlr --itsv --omd cat

# divided by the number of subpop strains
cat group.num.tsv |
    datamash transpose |
    sed 's/^Mut/group/' |
    tsv-select -H -f group,All |
    tsv-join -H -k group -f ../isolates/strain.tsv -a num |
    sed 1d |
    perl -nla -e '$num = $F[1]/$F[2];printf qq{%s\t%.3f\n}, $F[0], $num;' |
    sed '1igroup\tsnp_per_strain' \
    > group.perstrain.tsv

# plot
Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    data <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(data, aes(x = reorder(group, -snp_per_strain), y = snp_per_strain)) +
         geom_bar(stat="identity") +
         geom_text(aes(label = snp_per_strain), vjust=1.6, color="white", size=3.5)+
         theme(axis.text.x = element_text(angle = 315))
    ggsave(p, height = 6, width = 15, file = "../results/group.num.pdf")
' group.perstrain.tsv

# mv all grouped info to one file
rm group.snp.tsv

for group in $(cat ../isolates/group.lst)
do
    echo "==> ${group}"
    cat ${group}.snp.tsv |
        tsv-join -k 1,2,3,4 -f random.snp.tsv -a 5,6,7 |
        awk -v group=${group} '{print ($0 "\t" group)}' \
    >> group.snp.tsv
done

cat group.snp.tsv | wc -l
#893
# totally 893 snps (an snp could exist among multiple groups)

cat group.snp.tsv | tsv-filter --ge 5:0.05 | wc -l
#453
# 453 snps population freq >= 0.05

# uniq all snps
cat group.snp.tsv |
    tsv-select -f 1,2,3,4 |
    tsv-uniq |
    wc -l
#282
# 282 snps among groups

cat group.snp.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-select -f 1,2,3,4 |
    tsv-uniq |
    wc -l
#127
# totally 127 high freq (>= 0.05) snps found among groups

# extract for plot
cat group.snp.tsv |
    tsv-select -f 9,10,8,5,11 |
    perl -nla -e '
        print join("\t",@F) if $F[0] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[0] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1itype\tgene\tfit\tfreq\tgroup' \
    > ../results/group.tsv

mkdir ~/data/yeast/vcf/snp

# gene count
cat group.snp.tsv |
    tsv-summarize -g 1,2,3,4,10 --count |
    tsv-sort -r -nk 6,6 |
    tsv-filter --ge 6:10 |
    sed '1ichr\tpos\tREF\tALT\tgene\tnum' |
    mlr --itsv --omd cat

# freq among groups
cat group.snp.tsv |
    tsv-summarize -g 1,2,3,4,10 --count |
    tsv-sort -r -nk 6,6 |
    tsv-filter --ge 6:10 |
    parallel --colsep '\t' -j 4 -k '
        cat group.snp.tsv |
        tsv-filter --str-eq 1:{1} --eq 2:{2} --str-eq 3:{3} --str-eq 4:{4} |
        tsv-select -f 11,9,5 |
        sed "1igroup\ttype\tfreq" \
        > snp/{1}_{2}_{3}_{4}.tsv
    '

# plot
for file in $(ls snp/*.tsv)
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

SNPs among groups:

| Mut        | Bakery | Beer | Bioethanol | Cider | Clinical | Dairy | Distillery | Fermentation | Flower | Fruit | Human | Industrial | Insect | Lab_strain | Nature | Palm_wine | Probiotic | Sake | Soil | Tree | Unknown | Water | Wine |
|------------|--------|------|------------|-------|----------|-------|------------|--------------|--------|-------|-------|------------|--------|------------|--------|-----------|-----------|------|------|------|---------|-------|------|
| All        | 35     | 70   | 24         | 21    | 62       | 21    | 43         | 40           | 22     | 64    | 35    | 35         | 45     | 13         | 91     | 51        | 7         | 20   | 40   | 70   | 47      | 39    | 77   |
| Pos        | 35     | 70   | 24         | 21    | 62       | 21    | 43         | 40           | 22     | 63    | 35    | 35         | 45     | 13         | 91     | 51        | 7         | 20   | 40   | 70   | 47      | 39    | 77   |
| One_SNP    | 35     | 70   | 24         | 21    | 62       | 21    | 43         | 40           | 22     | 62    | 35    | 35         | 45     | 13         | 91     | 51        | 7         | 20   | 40   | 70   | 47      | 39    | 77   |
| Two_SNPs   | 0      | 0    | 0          | 0     | 0        | 0     | 0          | 0            | 0      | 1     | 0     | 0          | 0      | 0          | 0      | 0         | 0         | 0    | 0    | 0    | 0       | 0     | 0    |
| Three_SNPs | 0      | 0    | 0          | 0     | 0        | 0     | 0          | 0            | 0      | 0     | 0     | 0          | 0      | 0          | 0      | 0         | 0         | 0    | 0    | 0    | 0       | 0     | 0    |

SNPs gene count:

| chr          | pos     | REF | ALT | gene  | num |
|--------------|---------|-----|-----|-------|-----|
| chromosome12 | 609213  | A   | G   | EST1  | 23  |
| chromosome7  | 457631  | T   | C   | GET1  | 22  |
| chromosome7  | 457592  | A   | G   | GET1  | 22  |
| chromosome12 | 1006428 | T   | C   | TSR2  | 22  |
| chromosome6  | 223291  | G   | A   | RPL29 | 21  |
| chromosome5  | 70005   | A   | G   | IES6  | 21  |
| chromosome14 | 477742  | G   | A   | EOS1  | 21  |
| chromosome14 | 477705  | C   | T   | EOS1  | 21  |
| chromosome3  | 211204  | A   | G   | BUD23 | 20  |
| chromosome5  | 69939   | A   | G   | IES6  | 19  |
| chromosome12 | 369873  | T   | C   | CCW12 | 19  |
| chromosome12 | 609230  | C   | T   | EST1  | 18  |
| chromosome3  | 211292  | G   | C   | BUD23 | 16  |
| chromosome7  | 457619  | G   | A   | GET1  | 13  |
| chromosome6  | 223290  | T   | A   | RPL29 | 12  |
| chromosome5  | 70029   | T   | G   | IES6  | 11  |
| chromosome13 | 500539  | G   | C   | ASC1  | 11  |
| chromosome12 | 609150  | A   | G   | EST1  | 11  |
| chromosome7  | 394084  | C   | T   | RAD6  | 10  |
| chromosome5  | 69963   | C   | G   | IES6  | 10  |
| chromosome12 | 1006476 | G   | A   | TSR2  | 10  |

- Mean and median of fitness

```bash
cd ~/data/yeast/vcf

# fitness mean and median among groups
# plot
Rscript -e '
    library(ggplot2)
    library(readr)
    library(plyr)
    args <- commandArgs(T)
    data <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(data, aes(x = fit, fill = type)) +
         geom_histogram(alpha = 0.5, position = "identity") +
         facet_wrap(~group)
    ggsave(p, height = 6, width = 12, file = "../results/group.fit.pdf")
' ../results/group.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    data <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(data, aes(x = type, y = fit, fill = type))+
         geom_boxplot() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    ggsave(p, height = 6, width = 15, file = "../results/group.fit.box.pdf")
' ../results/group.tsv

# count
cat group.snp.tsv |
    tsv-summarize -g 9,11 --mean 8 --median 8 --count |
    sed '1itype\tgroup\tfit_mean\tfit_median\tcount' |
    mlr --itsv --omd cat
```

| type                   | group        | fit_mean       | fit_median     | count |
|------------------------|--------------|----------------|----------------|-------|
| Synonymous_mutation    | Bakery       | 0.984752522254 | 0.98724477975  | 16    |
| Nonsynonymous_mutation | Bakery       | 0.987965685909 | 0.98797192625  | 15    |
| Synonymous_mutation    | Beer         | 0.984413066858 | 0.98658440575  | 35    |
| Nonsynonymous_mutation | Beer         | 0.989262805563 | 0.9917333615   | 30    |
| Synonymous_mutation    | Bioethanol   | 0.986866381265 | 0.993277023    | 13    |
| Nonsynonymous_mutation | Bioethanol   | 0.986384899024 | 0.987598707    | 9     |
| Nonsynonymous_mutation | Cider        | 0.985603878747 | 0.98747437925  | 9     |
| Synonymous_mutation    | Cider        | 0.984326260382 | 0.98386353025  | 11    |
| Nonsynonymous_mutation | Clinical     | 0.988420517561 | 0.9884840435   | 24    |
| Synonymous_mutation    | Clinical     | 0.985689947094 | 0.987584157281 | 34    |
| Nonsynonymous_mutation | Dairy        | 0.985667778537 | 0.986215389375 | 6     |
| Synonymous_mutation    | Dairy        | 0.99061800815  | 0.9940471665   | 13    |
| Nonsynonymous_mutation | Distillery   | 0.990814705196 | 0.99272248975  | 14    |
| Synonymous_mutation    | Distillery   | 0.986273908723 | 0.98744413375  | 26    |
| Synonymous_mutation    | Fermentation | 0.984257335963 | 0.986211182614 | 24    |
| Nonsynonymous_mutation | Fermentation | 0.987159823702 | 0.987689188267 | 13    |
| Synonymous_mutation    | Flower       | 0.984196724407 | 0.987483555625 | 14    |
| Nonsynonymous_mutation | Flower       | 0.987912223069 | 0.987560714706 | 6     |
| Synonymous_mutation    | Fruit        | 0.987341498499 | 0.98857693925  | 40    |
| Nonsynonymous_mutation | Fruit        | 0.989193583458 | 0.987521812    | 20    |
| Nonsynonymous_mutation | Human        | 0.985104257254 | 0.987109888837 | 15    |
| Synonymous_mutation    | Human        | 0.984297448669 | 0.982130373375 | 16    |
| Nonsynonymous_mutation | Industrial   | 0.987222439578 | 0.987791252625 | 14    |
| Synonymous_mutation    | Industrial   | 0.986468714323 | 0.98807367875  | 19    |
| Nonsynonymous_mutation | Insect       | 0.991944864397 | 0.989291307    | 10    |
| Synonymous_mutation    | Insect       | 0.985738207509 | 0.9889152495   | 30    |
| Synonymous_mutation    | Lab_strain   | 0.981781640668 | 0.987104177062 | 9     |
| Nonsynonymous_mutation | Lab_strain   | 0.9838339215   | 0.9838339215   | 2     |
| Synonymous_mutation    | Nature       | 0.98789692653  | 0.98955296025  | 58    |
| Nonsynonymous_mutation | Nature       | 0.989460055308 | 0.990928906    | 29    |
| Synonymous_mutation    | Palm_wine    | 0.98571137937  | 0.98795345925  | 31    |
| Nonsynonymous_mutation | Palm_wine    | 0.987947942565 | 0.987598707    | 15    |
| Nonsynonymous_mutation | Probiotic    | 0.992545284917 | 0.9943413815   | 3     |
| Synonymous_mutation    | Probiotic    | 0.973690638649 | 0.97665198325  | 3     |
| Synonymous_mutation    | Sake         | 0.983589069915 | 0.984590859208 | 14    |
| Nonsynonymous_mutation | Sake         | 0.984269561988 | 0.9809705645   | 3     |
| Synonymous_mutation    | Soil         | 0.987892931407 | 0.991103581    | 21    |
| Nonsynonymous_mutation | Soil         | 0.991429979907 | 0.99420745     | 15    |
| Synonymous_mutation    | Tree         | 0.98577578033  | 0.9880641375   | 45    |
| Nonsynonymous_mutation | Tree         | 0.987213133054 | 0.987303140125 | 18    |
| Nonsynonymous_mutation | Unknown      | 0.988712827819 | 0.98797192625  | 19    |
| Synonymous_mutation    | Unknown      | 0.985133560791 | 0.9877505665   | 23    |
| Nonsynonymous_mutation | Water        | 0.98799604383  | 0.988614516    | 18    |
| Synonymous_mutation    | Water        | 0.98662720023  | 0.988260569875 | 18    |
| Synonymous_mutation    | Wine         | 0.987165814852 | 0.98913000425  | 38    |
| Nonsynonymous_mutation | Wine         | 0.98840118878  | 0.99077134075  | 35    |

- Mean and median of frequencies

```bash
cd ~/data/yeast/vcf

# plot
Rscript -e '
    library(ggplot2)
    library(readr)
    library(plyr)
    args <- commandArgs(T)
    data <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(data, aes(x = freq, fill = type)) +
         geom_histogram(alpha = 0.5, position = "identity") +
         facet_wrap(~group)
    ggsave(p, height = 6, width = 12, file = "../results/group.freq.pdf")
' ../results/group.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    data <- read_tsv(args[1], show_col_types = FALSE)
    p <- ggplot(data, aes(x = type, y = freq, fill = type))+
         geom_boxplot() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    ggsave(p, height = 6, width = 15, file = "../results/group.freq.box.pdf")
' ../results/group.tsv

# count
cat group.snp.tsv |
    tsv-summarize -g 9,11 --mean 5 --median 5 --count |
    sed '1itype\tgroup\tfreq_mean\tfreq_median\tcount' |
    mlr --itsv --omd cat
```

| type                   | group        | freq_mean       | freq_median | count |
|------------------------|--------------|-----------------|-------------|-------|
| Synonymous_mutation    | Bakery       | 0.29983110625   | 0.2027025   | 16    |
| Nonsynonymous_mutation | Bakery       | 0.09770914      | 0.0675676   | 15    |
| Synonymous_mutation    | Beer         | 0.166705979429  | 0.0423729   | 35    |
| Nonsynonymous_mutation | Beer         | 0.0589567486667 | 0.0338983   | 30    |
| Synonymous_mutation    | Bioethanol   | 0.4059827       | 0.185185    | 13    |
| Nonsynonymous_mutation | Bioethanol   | 0.148148177778  | 0.0555556   | 9     |
| Nonsynonymous_mutation | Cider        | 0.216298922222  | 0.205882    | 9     |
| Synonymous_mutation    | Cider        | 0.361463854545  | 0.117647    | 11    |
| Nonsynonymous_mutation | Clinical     | 0.0712669529167 | 0.0211162   | 24    |
| Synonymous_mutation    | Clinical     | 0.136218309118  | 0.02816965  | 34    |
| Nonsynonymous_mutation | Dairy        | 0.370370266667  | 0.1944443   | 6     |
| Synonymous_mutation    | Dairy        | 0.410826238462  | 0.111111    | 13    |
| Nonsynonymous_mutation | Distillery   | 0.113300485714  | 0.0344828   | 14    |
| Synonymous_mutation    | Distillery   | 0.239342519231  | 0.0603448   | 26    |
| Synonymous_mutation    | Fermentation | 0.307291675     | 0.1875001   | 24    |
| Nonsynonymous_mutation | Fermentation | 0.146367523077  | 0.0277778   | 13    |
| Synonymous_mutation    | Flower       | 0.3571429       | 0.214286    | 14    |
| Nonsynonymous_mutation | Flower       | 0.202380916667  | 0.0892858   | 6     |
| Synonymous_mutation    | Fruit        | 0.1377848175    | 0.0425532   | 40    |
| Nonsynonymous_mutation | Fruit        | 0.065957455     | 0.02659575  | 20    |
| Nonsynonymous_mutation | Human        | 0.209677433333  | 0.0322581   | 15    |
| Synonymous_mutation    | Human        | 0.4153226375    | 0.177419    | 16    |
| Nonsynonymous_mutation | Industrial   | 0.0955782714286 | 0.06904765  | 14    |
| Synonymous_mutation    | Industrial   | 0.2616541       | 0.183333    | 19    |
| Nonsynonymous_mutation | Insect       | 0.11            | 0.05        | 10    |
| Synonymous_mutation    | Insect       | 0.2225          | 0.125       | 30    |
| Synonymous_mutation    | Lab_strain   | 0.611111111111  | 0.5         | 9     |
| Nonsynonymous_mutation | Lab_strain   | 0.5             | 0.5         | 2     |
| Synonymous_mutation    | Nature       | 0.0921880668966 | 0.0192308   | 58    |
| Nonsynonymous_mutation | Nature       | 0.0514511068966 | 0.0192308   | 29    |
| Synonymous_mutation    | Palm_wine    | 0.258064529032  | 0.0666667   | 31    |
| Nonsynonymous_mutation | Palm_wine    | 0.12222224      | 0.0333333   | 15    |
| Nonsynonymous_mutation | Probiotic    | 0.5             | 0.5         | 3     |
| Synonymous_mutation    | Probiotic    | 1               | 1           | 3     |
| Synonymous_mutation    | Sake         | 0.6398177       | 0.9787235   | 14    |
| Nonsynonymous_mutation | Sake         | 0.386524866667  | 0.148936    | 3     |
| Synonymous_mutation    | Soil         | 0.247493747619  | 0.0921053   | 21    |
| Nonsynonymous_mutation | Soil         | 0.10087722      | 0.0789474   | 15    |
| Synonymous_mutation    | Tree         | 0.143515788889  | 0.03125     | 45    |
| Nonsynonymous_mutation | Tree         | 0.0744391944444 | 0.03125     | 18    |
| Nonsynonymous_mutation | Unknown      | 0.0817669210526 | 0.0357143   | 19    |
| Synonymous_mutation    | Unknown      | 0.221273278261  | 0.0357143   | 23    |
| Nonsynonymous_mutation | Water        | 0.0862572888889 | 0.0526316   | 18    |
| Synonymous_mutation    | Water        | 0.269736877778  | 0.0789473   | 18    |
| Synonymous_mutation    | Wine         | 0.0897177571053 | 0.00806452  | 38    |
| Nonsynonymous_mutation | Wine         | 0.0366996411429 | 0.00403226  | 35    |

- Linear analysis

```bash
cd ~/data/yeast/vcf

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
data <- read_tsv(args[1], show_col_types = FALSE)
p <- ggplot(data, aes(x = freq, y = fit, color = group)) +
     geom_point() +
     geom_smooth(method = "lm") +
     facet_wrap(~group)
ggsave(p, height = 6, width = 15, file = "../results/group.lm.pdf")
' ../results/group.tsv
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

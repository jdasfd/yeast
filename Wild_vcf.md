# Check wild yeast vcf

The [article](https://www.nature.com/articles/s41586-022-04823-w) contains thousands of mutations from 21 genes (synonymous or non-synonymous mutations). The goal is to validate all mutations among wild yeasts.

`mkdir -p ~/data/yeast/info`

Download all supplementary files and tables from the [website](https://www.nature.com/articles/s41586-022-04823-w) and save them into info dir.

Steps followed here are processes to convert the mutation info into vcf format for better analyze.

There are some data from the [README.md](https://github.com/wang-q/pars#readme). Better completing those steps.

Scripts could be partly found in the github repo [pars](https://github.com/wang-q/pars). Make sure you clone it before.

## Data preparation

### Split strains from 1002 genomes project according to ecological groups

- Acquire info and split to subpopulations

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

- Split the `1011Matrix.gvcf.gz` into small files according to subpopulations

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
```

### Extract 21 gene names and sequences

According to SGD, there is system name and standard name for one same gene. The standard name was adopted in the article. When it comes to different processes, two names were respectively useful in different aspects.

After that, the sequences were extracted from the mutation info.

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
        perl -nae 'print $F[2]; END{print qq{\n};}' \
        >> gene.mut.fa
done

# All genes and their length
faops size gene.mut.fa
#ADA2    150
#ASC1    142
#BFR1    150
#BUD23   150
#CCW12   150
#EOS1    150
#EST1    150
#GET1    150
#GIM5    150
#IES6    150
#LSM1    112
#PAF1    150
#PRS3    149
#RAD6    149
#RPL29   150
#RPL39   146
#RPS7A   141
#SNF6    150
#TSR2    150
#VMA21   147
#VMA7    119
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
perl ~/Scripts/pars/blastn_transcript.pl -i 90 -f gene.blast -m 0

cat gene.blast.tsv | wc -l
#19
# LSM1 and VMA7 excluded
# too short (showed in previous step)
```

### Get vcf of genes and of subpopulations

There are few things should be considered:

1. The file format of 1002 genome project VCF is [VCF v4.1](http://samtools.github.io/hts-specs/VCFv4.1.pdf). In this format, all  info of mutations was recorded on the forward strand of reference yeast genome R64.

2. In the article, each position of mutations was recorded according to the coding region. So after I researched on the [genome browser](https://browse.yeastgenome.org/), it should be noticed that the blast result should be relocated on the forward strand.

In other words, for the forward strand (+), the genome location is alright, but it should be upside down for the reverse strand (-).

- Forward strand: left 5' -> right 3'.
- Reverse strand: right 5' -> left 3'.

```bash
cd ~/data/yeast

cat info/fit.tsv |
    tsv-filter --str-ne 1:LSM1 --str-ne 1:VMA7 \
    > gene/fit.filter.tsv

perl scripts/loc2vcf.pl -b gene/gene.blast.tsv \
    -t gene/fit.filter.tsv |
    tsv-sort -nk 2,2 \
    > vcf/gene.vcf

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
        > vcf/region.1based.bed

# combine all vcf and info into 1 file
for group in $(cat isolates/group.lst)
do
    bcftools view -Ov -R vcf/region.1based.bed \
                  -v snps,mnps,ref --threads 10 \
                  vcf/group/1011Matrix.${group}.bcf |
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF{1}\t%AC{1}\t%AN{1}\n' |
    tsv-join -f vcf/gene.vcf -k 1,2,3,4 -a 5,6,7 -z |
    tsv-filter --ne 5:0 |
    awk -v col="${group}" '{print ($0 "\t" col)}' \
    >> vcf/all.vcf.tsv
done

# all mutations undetected
perl scripts/loc2vcf.pl -b gene/gene.blast.tsv \
    -t gene/fit.filter.tsv |
    tsv-join -f vcf/all.vcf.tsv -k 1,2,3,4 -e \
    > vcf/other.vcf.tsv

cat vcf/other.vcf.tsv | wc -l
#7756 (8004 - 248)
```

## Statistical analysis

### Basic info

- The numbers of SNPs occurred among subpopulations

```bash
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
mv tmp group_num.tsv

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
    ggsave(p, height = 6, width = 15, file = "group_num.pdf")
' group_num.tsv
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

cat all.vcf.tsv |
    tsv-summarize -g 1,2,3,4,10 --count |
    wc -l
#248

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
# output into a tsv for plot
cat all.vcf.tsv |
    tsv-summarize -g 11,9 --count --mean 5 --median 5 |
    perl -nla -e '
        print join("\t",@F) if $F[1] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[1] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1igroup\ttype\tnum\tfreq_mean\tfreq_median' \
    > group_freq.tsv

cat group_freq.tsv | mlr --itsv --omd cat

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
    ggsave(p1, height = 6, width = 15, file = "group_freq_mean.pdf")
    ggsave(p2, height = 6, width = 15, file = "group_freq_median.pdf")
' group_freq.tsv

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-summarize -g 11,9 --count --mean 5 --median 5 |
    perl -nla -e '
        print join("\t",@F) if $F[1] =~ s/^Non.+$/N_mut/;
        print join("\t",@F) if $F[1] =~ s/^Sy.+$/S_mut/;
    ' |
    sed '1igroup\ttype\tnum\tfreq_mean\tfreq_median' \
    > group_high_freq.tsv

cat group_high_freq.tsv | mlr --itsv --omd cat

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
    ggsave(p1, height = 6, width = 15, file = "group_high_freq_mean.pdf")
    ggsave(p2, height = 6, width = 15, file = "group_high_freq_median.pdf")
' group_high_freq.tsv
```

| group        | type  | num | freq_mean       | freq_median |
|--------------|-------|-----|-----------------|-------------|
| Bakery       | S_mut | 15  | 0.31441444      | 0.27027     |
| Bakery       | N_mut | 14  | 0.102757864286  | 0.0675676   |
| Beer         | S_mut | 33  | 0.175525345455  | 0.0423729   |
| Beer         | N_mut | 26  | 0.0634637838462 | 0.0338983   |
| Bioethanol   | S_mut | 12  | 0.430555341667  | 0.304131    |
| Bioethanol   | N_mut | 8   | 0.136574075     | 0.03703705  |
| Cider        | N_mut | 9   | 0.216298922222  | 0.205882    |
| Cider        | S_mut | 10  | 0.38823524      | 0.2058825   |
| Clinical     | N_mut | 23  | 0.0739591773913 | 0.0233645   |
| Clinical     | S_mut | 31  | 0.146980374839  | 0.0280374   |
| Dairy        | N_mut | 6   | 0.370370266667  | 0.1944443   |
| Dairy        | S_mut | 13  | 0.410826238462  | 0.111111    |
| Distillery   | N_mut | 13  | 0.119363384615  | 0.0344828   |
| Distillery   | S_mut | 25  | 0.2461576       | 0.0517241   |
| Fermentation | S_mut | 22  | 0.332702027273  | 0.277778    |
| Fermentation | N_mut | 12  | 0.15625         | 0.0277778   |
| Flower       | S_mut | 13  | 0.368131892308  | 0.214286    |
| Flower       | N_mut | 6   | 0.202380916667  | 0.0892858   |
| Fruit        | S_mut | 36  | 0.150369286111  | 0.0425532   |
| Fruit        | N_mut | 18  | 0.0721040277778 | 0.0319149   |
| Human        | N_mut | 14  | 0.221198171429  | 0.02419355  |
| Human        | S_mut | 15  | 0.43870974      | 0.193548    |
| Industrial   | N_mut | 14  | 0.0955782714286 | 0.06904765  |
| Industrial   | S_mut | 18  | 0.274338588889  | 0.2         |
| Insect       | N_mut | 9   | 0.116666666667  | 0.05        |
| Insect       | S_mut | 28  | 0.227678571429  | 0.1         |
| Lab_strain   | S_mut | 8   | 0.625           | 0.5         |
| Lab_strain   | N_mut | 1   | 0.5             | 0.5         |
| Nature       | S_mut | 53  | 0.0990637173585 | 0.0192308   |
| Nature       | N_mut | 30  | 0.0497360686667 | 0.0192308   |
| Palm_wine    | S_mut | 26  | 0.290384630769  | 0.0666667   |
| Palm_wine    | N_mut | 11  | 0.148484872727  | 0.0333333   |
| Probiotic    | N_mut | 3   | 0.5             | 0.5         |
| Probiotic    | S_mut | 3   | 1               | 1           |
| Sake         | S_mut | 13  | 0.688216115385  | 0.989362    |
| Sake         | N_mut | 3   | 0.386524866667  | 0.148936    |
| Soil         | S_mut | 20  | 0.255921065     | 0.13815815  |
| Soil         | N_mut | 13  | 0.106275330769  | 0.0789474   |
| Tree         | S_mut | 41  | 0.153325256098  | 0.03125     |
| Tree         | N_mut | 16  | 0.08081440625   | 0.03125     |
| Unknown      | N_mut | 19  | 0.0808270684211 | 0.0357143   |
| Unknown      | S_mut | 20  | 0.251785705     | 0.0714286   |
| Water        | N_mut | 15  | 0.09298244      | 0.0526316   |
| Water        | S_mut | 18  | 0.269736877778  | 0.0789473   |
| Wine         | S_mut | 31  | 0.108155583871  | 0.00806452  |
| Wine         | N_mut | 32  | 0.0402662403125 | 0.007056455 |

| group        | type  | num | freq_mean      | freq_median |
|--------------|-------|-----|----------------|-------------|
| Bakery       | S_mut | 11  | 0.420147463636 | 0.364865    |
| Bakery       | N_mut | 11  | 0.127097236364 | 0.0810811   |
| Beer         | N_mut | 9   | 0.1447334      | 0.059322    |
| Beer         | S_mut | 15  | 0.357907153333 | 0.245763    |
| Bioethanol   | N_mut | 4   | 0.25462965     | 0.2592595   |
| Bioethanol   | S_mut | 8   | 0.6319441375   | 0.7129625   |
| Cider        | N_mut | 9   | 0.216298922222 | 0.205882    |
| Cider        | S_mut | 8   | 0.4779411      | 0.3088235   |
| Clinical     | N_mut | 7   | 0.202879742857 | 0.122642    |
| Clinical     | S_mut | 13  | 0.333597338462 | 0.264423    |
| Dairy        | N_mut | 4   | 0.5370369      | 0.583333    |
| Dairy        | S_mut | 9   | 0.576954788889 | 0.740741    |
| Distillery   | N_mut | 5   | 0.2586206      | 0.241379    |
| Distillery   | S_mut | 13  | 0.446854769231 | 0.517241    |
| Fermentation | S_mut | 15  | 0.474074073333 | 0.472222    |
| Fermentation | N_mut | 5   | 0.34166664     | 0.25        |
| Flower       | S_mut | 13  | 0.368131892308 | 0.214286    |
| Flower       | N_mut | 5   | 0.23571424     | 0.107143    |
| Fruit        | N_mut | 5   | 0.19148938     | 0.106383    |
| Fruit        | S_mut | 16  | 0.30563418125  | 0.209528    |
| Human        | N_mut | 6   | 0.49193555     | 0.564516    |
| Human        | S_mut | 10  | 0.6467743      | 0.806452    |
| Industrial   | N_mut | 9   | 0.130158811111 | 0.1         |
| Industrial   | S_mut | 14  | 0.344387714286 | 0.2416665   |
| Insect       | N_mut | 7   | 0.142857142857 | 0.05        |
| Insect       | S_mut | 26  | 0.243269230769 | 0.125       |
| Lab_strain   | S_mut | 8   | 0.625          | 0.5         |
| Lab_strain   | N_mut | 1   | 0.5            | 0.5         |
| Nature       | N_mut | 7   | 0.145846857143 | 0.0882353   |
| Nature       | S_mut | 16  | 0.28485573125  | 0.216346    |
| Palm_wine    | N_mut | 5   | 0.29000008     | 0.166667    |
| Palm_wine    | S_mut | 17  | 0.426470629412 | 0.183333    |
| Probiotic    | N_mut | 3   | 0.5            | 0.5         |
| Probiotic    | S_mut | 3   | 1              | 1           |
| Sake         | S_mut | 9   | 0.988179777778 | 1           |
| Sake         | N_mut | 2   | 0.569149       | 0.569149    |
| Soil         | S_mut | 19  | 0.268005552632 | 0.184211    |
| Soil         | N_mut | 8   | 0.157894775    | 0.0789474   |
| Tree         | S_mut | 18  | 0.318424888889 | 0.230469    |
| Tree         | N_mut | 4   | 0.250992       | 0.09375     |
| Unknown      | N_mut | 6   | 0.1904762      | 0.142857    |
| Unknown      | S_mut | 11  | 0.435064927273 | 0.428571    |
| Water        | S_mut | 16  | 0.3001645125   | 0.131579    |
| Water        | N_mut | 11  | 0.117224854545 | 0.0789474   |
| Wine         | N_mut | 6   | 0.169054333333 | 0.09173405  |
| Wine         | S_mut | 4   | 0.75201625     | 0.947581    |

- The mean and median of all mutation fitness

```bash
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
cd ~/data/yeast/fitness

cat all.fit.tsv |
    tsv-filter --str-ne 7:Nonsense_mutation |
    tsv-select -f 1,6,7 |
    perl -nla -e '
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Non.+$/N_mut/;
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Sy.+$/S_mut/;
        ' |
    sed '1igroup\tfit\ttype' \
    > fitness.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1])
p <- ggplot(fit, aes(x = type, y = fit)) +
     geom_boxplot() +
     facet_grid(~group) +
     ylim(0.96, 1.01) +
     theme(axis.text.x = element_text(angle = 315))
ggsave(p, height = 6, width = 15, file = "fitness.pdf")
' fitness.tsv

cat all.fit.tsv |
    tsv-filter --str-ge 8:0.05 --str-ne 7:Nonsense_mutation |
    tsv-select -f 1,6,7 \
    > fitness_5.tsv

cat fitness_5.tsv \
<(cat all.fit.tsv |
      tsv-filter --str-ne 7:Nonsense_mutation |
      tsv-filter --str-eq 1:other |
      tsv-select -f 1,6,7) |
    perl -nla -e '
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Non.+$/N_mut/;
        print qq($F[0]\t$F[1]\t$F[2]) if $F[2] =~ s/^Sy.+$/S_mut/;
        ' |
    sed '1igroup\tfit\ttype' \
      > tmp && mv tmp fitness_5.tsv 

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1])
p <- ggplot(fit, aes(x = type, y = fit)) +
     geom_boxplot() +
     facet_grid(~group) +
     ylim(0.96, 1.01) +
     theme(axis.text.x = element_text(angle = 315))
ggsave(p, height = 6, width = 15, file = "fitness_5.pdf")
' fitness_5.tsv

# distribution of all muts fitness
cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-uniq -f 2,3,4,5 |
    tsv-select -f 6,7 |
    sed '1ifit\ttype' \
    > dist.tsv

Rscript -e '
    library(ggplot2)
    library(readr)
    args <- commandArgs(T)
    fit <- read_tsv(args[1])
    p <- ggplot(fit, aes(x = fit)) +
         geom_histogram() +
         facet_grid(~type)
    ggsave(p, height = 6, width = 12, file = "dist.pdf")
' dist.tsv
```

### Chi-square

The experiment from the original article was almost a simulation of random mutations. All detected mutations were not biased to either nonsynonymous mutation or synonymous mutation. So all unique fixed SNPs after would be random as well if there were no other reasons.

```bash
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

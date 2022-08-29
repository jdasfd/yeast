# Check wild yeast vcf

The [article](https://www.nature.com/articles/s41586-022-04823-w) contains thousands of mutations from 21 genes (synonymous or non-synonymous mutations). The goal is to validate all mutations among wild yeasts.

`mkdir -p ~/data/yeast/info`

Download all supplementary files and tables from the [website](https://www.nature.com/articles/s41586-022-04823-w) and save them into info dir.

Steps followed here are processes to convert the mutation info into vcf format for better analyze.

There are some data from the [README.md](https://github.com/wang-q/pars#readme). Better completing those steps.

Scripts could be partly found in the github repo [pars](https://github.com/wang-q/pars). Make sure you clone it before.

## Data preparation

### Split strains from 1002 genomes project according to ecological groups

- Acquire info and split to subgroups

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

- Split the `1011Matrix.gvcf.gz` into small files according to subgroups

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

### Get vcf of genes and subpopulations

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
```

## Statistical analysis

### Basic info

```bash
cd ~/data/yeast/vcf

cat all.vcf.tsv | wc -l
#819
# totally 819 snps

cat all.vcf.tsv | tsv-filter --ge 5:0.05 | wc -l
#431
# 431 snps population freq >= 0.05

cat all.vcf.tsv |
    tsv-sort -k 1,1 -k 2,2 -k 3,3 -k 4,4 |
    tsv-select -f 1,2,3,4 |
    tsv-uniq |
    wc -l
#248
# totally 248 snps found among wild groups

cat all.vcf.tsv |
    tsv-filter --ge 5:0.05 |
    tsv-sort -k 1,1 -k 2,2 -k 3,3 -k 4,4 |
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
    data <- read_tsv(args[1])
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

```bash
cd ~/data/yeast/vcf

cat ../isolates/group.lst | wc -l
#23

cat all.vcf.tsv |
    tsv-summarize -g 1,2,3,4,10 --count |
    wc -l
#248

# all snps repeated more than 10 times
cat all.vcf.tsv |
    tsv-summarize -g 1,2,3,4,10 --count |
    tsv-sort -r -nk 6,6 |
    tsv-select -f 1,5,6 |
    tsv-filter --ge 3:10 |
    sed '1ichr\tgene\tgroup' |
    mlr --itsv --omd cat
```

| chr | gene | group |
| --- | --- | --- |
| chromosome12 | EST1 | 23 |
| chromosome7 | GET1 | 22 |
| chromosome7 | GET1 | 22 |
| chromosome12 | TSR2 | 22 |
| chromosome6 | RPL29 | 21 |
| chromosome5 | IES6 | 21 |
| chromosome14 | EOS1 | 21 |
| chromosome14 | EOS1 | 21 |
| chromosome3 | BUD23 | 20 |
| chromosome5 | IES6 | 19 |
| chromosome12 | CCW12 | 19 |
| chromosome12 | EST1 | 18 |
| chromosome3 | BUD23 | 16 |
| chromosome7 | GET1 | 13 |
| chromosome6 | RPL29 | 12 |
| chromosome5 | IES6 | 11 |
| chromosome13 | ASC1 | 11 |
| chromosome12 | EST1 | 11 |
| chromosome5 | IES6 | 10 |
| chromosome12 | TSR2 | 10 |

```bash
cd ~/data/yeast/vcf

# freq among groups
cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-select -f 2,3,4,5 |
    tsv-uniq |
    parallel --colsep '\t' -j 4 -k '
        cat all.fit.tsv |
        tsv-filter --str-eq 2:{1} --eq 3:{2} --str-eq 4:{3} --str-eq 5:{4} |
        tsv-select -f 1,7,8 |
        sed "1igroup\ttype\tfreq" \
        > freq/{1}_{2}_{3}_{4}.tsv
    '

cd ~/data/yeast/fitness/freq

for name in $(wc -l *.tsv |
              grep -v 'total$' |
              datamash reverse -W |
              tsv-filter --eq 2:2 |
              tsv-select -f 1)
do
    rm ${name}
done

for file in $(ls *.tsv)
do
    echo "==>${file}"
    Rscript -e '
        library(ggplot2)
        library(readr)
        args <- commandArgs(T)
        freq <- read_tsv(args[1], show_col_types = FALSE)
        save <- paste0(args[1], ".pdf")
        p <- ggplot(freq, aes(x = group, y = freq))+
             geom_col() +
             theme(axis.text.x = element_text(angle = 315))
        ggsave(p, height = 6, width = 10, file = save)
    ' ${file}
done

cd ~/data/yeast/fitness
# freq mean and median
# snps existent in wild groups
cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-uniq -f 2,3,4,5 |
    tsv-summarize -g 7 --count --mean 8 --median 8 |
    sed '1itype\tnum\tfreq_mean\tfreq_median' |
    mlr --itsv --omd cat

cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-summarize -g 1,7 --count --mean 8 --median 8 |
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
    p1 <- ggplot(freq, aes(x = type, y = freq_mean))+
         geom_col() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    p2 <- ggplot(freq, aes(x = type, y = freq_median))+
         geom_col() +
         theme(axis.text.x = element_text(angle = 315)) +
         facet_grid(~group)
    ggsave(p1, height = 6, width = 15, file = "group_freq_mean.pdf")
    ggsave(p2, height = 6, width = 15, file = "group_freq_median.pdf")
' group_freq.tsv

# fitness mean and median
# uniq snps from different groups
# snps existent in wild groups
cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-uniq -f 2,3,4,5 |
    tsv-summarize -g 7 --count --mean 6 --median 6 |
    sed '1itype\tnum\tfit_mean\tfit_median' |
    mlr --itsv --omd cat

# snps nonexistent in wild groups
cat all.fit.tsv |
    tsv-filter --str-eq 1:other |
    tsv-summarize -g 7 --count --mean 6 --median 6 |
    sed '1itype\tnum\tfit_mean\tfit_median' |
    mlr --itsv --omd cat
```

repeated snps:

| TSR2  | XII  | 1006428 | T   | C   | 22  |
|-------|------|---------|-----|-----|-----|
| IES6  | V    | 69766   | A   | G   | 22  |
| BFR1  | XV   | 718752  | T   | C   | 21  |
| EST1  | XII  | 607457  | T   | C   | 21  |
| RPL39 | X    | 76457   | A   | T   | 5   |
| CCW12 | XII  | 369768  | T   | C   | 4   |
| BUD23 | III  | 210821  | G   | A   | 4   |
| RPS7A | XV   | 506434  | A   | G   | 3   |
| RPL39 | X    | 76435   | T   | C   | 3   |
| RPL39 | X    | 76417   | C   | T   | 3   |
| PRS3  | VIII | 80765   | C   | T   | 3   |
| SNF6  | VIII | 54942   | T   | A   | 3   |
| GET1  | VII  | 457187  | G   | A   | 3   |
| BUD23 | III  | 210807  | G   | A   | 3   |
| BFR1  | XV   | 718726  | T   | A   | 2   |
| TSR2  | XII  | 1006524 | C   | G   | 2   |
| TSR2  | XII  | 1006437 | T   | A   | 2   |
| LSM1  | X    | 187267  | G   | A   | 2   |
| VMA21 | VII  | 698670  | C   | T   | 2   |
| ADA2  | IV   | 1356180 | C   | T   | 2   |
| ADA2  | IV   | 1356133 | T   | C   | 2   |

existent snps fit:

| type                   | num | fit_mean       | fit_median     |
|------------------------|-----|----------------|----------------|
| Nonsynonymous_mutation | 33  | 0.987204450667 | 0.98809565525  |
| Synonymous_mutation    | 30  | 0.990295101505 | 0.99229843675  |
| Nonsense_mutation      | 3   | 0.92759421177  | 0.904949830064 |

nonexistent snps fit:

| type                   | num  | fit_mean       | fit_median     |
|------------------------|------|----------------|----------------|
| Nonsynonymous_mutation | 6273 | 0.984944235527 | 0.98805154375  |
| Synonymous_mutation    | 1836 | 0.987772083901 | 0.988734612306 |
| Nonsense_mutation      | 166  | 0.934085408798 | 0.939873592125 |

existent snps freq:

| type                   | num | freq_mean       | freq_median |
|------------------------|-----|-----------------|-------------|
| Nonsynonymous_mutation | 33  | 0.0385074172727 | 0.0212766   |
| Synonymous_mutation    | 30  | 0.0821522103333 | 0.0344828   |
| Nonsense_mutation      | 3   | 0.0585839666667 | 0.0789474   |

| group        | type  | num | freq_mean       | freq_median |
|--------------|-------|-----|-----------------|-------------|
| Bakery       | N_mut | 3   | 0.1846847       | 0.216216    |
| Bakery       | S_mut | 2   | 0.466216        | 0.466216    |
| Beer         | N_mut | 7   | 0.0932202971429 | 0.0508475   |
| Beer         | S_mut | 6   | 0.217513996667  | 0.08050825  |
| Bioethanol   | N_mut | 2   | 0.277778        | 0.277778    |
| Bioethanol   | S_mut | 4   | 0.259259125     | 0.0185185   |
| Cider        | N_mut | 2   | 0.1911765       | 0.1911765   |
| Cider        | S_mut | 2   | 0.2794115       | 0.2794115   |
| Clinical     | N_mut | 6   | 0.101584015     | 0.016355145 |
| Clinical     | S_mut | 4   | 0.236587375     | 0.1577544   |
| Dairy        | N_mut | 3   | 0.475308666667  | 0.388889    |
| Dairy        | S_mut | 2   | 0.5314815       | 0.5314815   |
| Distillery   | S_mut | 6   | 0.160919566667  | 0.1034484   |
| Distillery   | N_mut | 5   | 0.1896552       | 0.0344828   |
| Fermentation | S_mut | 6   | 0.240740616667  | 0.1666665   |
| Fermentation | N_mut | 3   | 0.4027776       | 0.333333    |
| Flower       | N_mut | 2   | 0.285714        | 0.285714    |
| Flower       | S_mut | 2   | 0.3571425       | 0.3571425   |
| Fruit        | N_mut | 4   | 0.16489355      | 0.0851063   |
| Fruit        | S_mut | 6   | 0.160471316667  | 0.04397165  |
| Human        | N_mut | 2   | 0.44354855      | 0.44354855  |
| Human        | S_mut | 5   | 0.43225804      | 0.403226    |
| Industrial   | N_mut | 3   | 0.183333433333  | 0.25        |
| Industrial   | S_mut | 3   | 0.3222221       | 0.1         |
| Insect       | N_mut | 3   | 0.308333333333  | 0.225       |
| Insect       | S_mut | 2   | 0.6375          | 0.6375      |
| Lab_strain   | N_mut | 2   | 0.5             | 0.5         |
| Lab_strain   | S_mut | 2   | 0.75            | 0.75        |
| Nature       | N_mut | 8   | 0.0721625275    | 0.01442309  |
| Nature       | S_mut | 6   | 0.139423056667  | 0.01923079  |
| Palm_wine    | N_mut | 4   | 0.4124999       | 0.40833315  |
| Palm_wine    | S_mut | 5   | 0.42333352      | 0.166667    |
| Sake         | N_mut | 4   | 0.7202244       | 0.9298105   |
| Sake         | S_mut | 2   | 0.994681        | 0.994681    |
| Soil         | N_mut | 3   | 0.179824433333  | 0.0789474   |
| Soil         | S_mut | 4   | 0.2302632       | 0.2039472   |
| Soil         | N_mut | 2   | 0.0789474       | 0.0789474   |
| Tree         | N_mut | 6   | 0.172433166667  | 0.0546875   |
| Tree         | S_mut | 6   | 0.1961805       | 0.0234375   |
| Unknown      | N_mut | 4   | 0.202711725     | 0.1375663   |
| Unknown      | N_mut | 1   | 0.0178571       | 0.0178571   |
| Unknown      | S_mut | 5   | 0.1607143       | 0.0357143   |
| Water        | N_mut | 2   | 0.1578945       | 0.1578945   |
| Water        | S_mut | 2   | 0.394737        | 0.394737    |
| Wine         | N_mut | 11  | 0.0137648890909 | 0.00403226  |
| Wine         | N_mut | 1   | 0.00201613      | 0.00201613  |
| Wine         | S_mut | 4   | 0.03931439      | 0.019153215 |

- Chi-square

The experiment from the original article was almost a simulation of random mutations. All detected mutations were not biased to either nonsynonymous mutation or synonymous mutation. So fixed SNPs after uniquified would be random as well if there were no other reasons.

```bash
raw1=$(cat all.fit.tsv |
          tsv-filter --str-ne 1:other |
          tsv-uniq -f 2,3,4,5 |
          tsv-summarize -g 7 --count |
          tsv-filter --str-ne 1:Nonsense_mutation |
          datamash transpose |
          sed 1d)

raw2=$(cat all.fit.tsv |
          tsv-filter --str-eq 1:other |
          tsv-summarize -g 7 --count |
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
#[1,]   33   30
#[2,] 6273 1836
#
#        Pearson's Chi-squared test with Yates' continuity correction
#
#data:  x
#X-squared = 20.74, df = 1, p-value = 5.262e-06
```

- Fitness from grouped SNPs

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

## Original article

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

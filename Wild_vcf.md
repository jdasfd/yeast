# Check wild yeast vcf

The [article](https://www.nature.com/articles/s41586-022-04823-w) contains thousands of mutations from 21 genes (synonymous or non-synonymous mutations). The goal is to validate all mutations among wild yeasts.

`mkdir -p ~/data/yeast/info`

Download all supplementary files and tables from the [website](https://www.nature.com/articles/s41586-022-04823-w) and save them into info dir.

Steps followed here are processes to convert the mutation info into vcf format for better analyze.

There are some data from the [README.md](https://github.com/wang-q/pars#readme). Better completing those steps.

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
# result is 1011 after subtraction, split correctly
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

# formatdb
makeblastdb -dbtype nucl -in S288c.fa -parse_seqids

# blast every transcripts against genome
blastn -task blastn -evalue 1e-3 -num_threads 4 -num_descriptions 10 -num_alignments 10 -outfmt 0 \
    -dust yes -soft_masking true \
    -db S288c.fa -query gene.mut.fa -out gene.blast

perl ~/Scripts/pars/blastn_transcript.pl -f gene.blast -m 0

# extract CDS region
# because of the coding region was selected from the SGD coding region
faops some ../mrna-structure/sgd/orf_genomic_all.fasta \
    gene/sysname.lst gene/21orf.fa

# rename by standard name
faops replace -l 0 gene/21orf.fa \
    <(cat gene/std_sysname.tsv | tsv-select -f 2,1) \
    gene/std_orf.fa

rm gene/21orf.fa

# align to the CDS
for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"
    
    cat  gene/${gene}/${gene}.mut.fa \
        <(faops one -l 0 gene/std_orf.fa ${gene} stdout) |
        muscle -out gene/${gene}/${gene}.aln.fa -quiet
done

<!--
```bash
# extract CDS genomic location
cat ../mrna-structure/sgd/saccharomyces_cerevisiae.gff |
    perl -nla -e '
    next if /^#/;
    next if /^[ATCG]/;
    if ($F[2] eq CDS){
        $F[8] =~ /Name=(.+)_CDS/;
        $name = $1;
        print "$F[0]\t$F[3]\t$F[4]\t$name";
        }
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
        -t gene/${gene}/${gene}.tsv \
        -a gene/${gene}/${gene}.aln.fa \
        > gene/${gene}/${gene}.mut.vcf
done
```
-->

### Filtering with bcf of different groups using region

```bash
mkdir ~/data/yeast/vcf/region
cd ~/data/yeast/vcf

# region containing gene name
rm region/region_name.tsv

for gene in $(cat ../gene/stdname.lst)
do
    cat ../gene/${gene}/${gene}_region.tsv >> region/region_name.tsv
done

cat region/region_name.tsv |
    tsv-join -k 1 --filter-file ../gene/sys_stdname.tsv \
    --append-fields 2 |
    tsv-select -f 5,1,2,3,4 > tmp \
    && mv tmp region/region_name.tsv

# .bed file for bcftools
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

rm ../vcf/region/all.mut.vcf

# transfer all mutation to 1 vcf
for gene in $(cat ../gene/stdname.lst)
do
    cat ../gene/${gene}/${gene}.mut.vcf \
        >> ../vcf/region/all.mut.vcf
done

cat region/all.mut.vcf | wc -l
#8341
# all mutation with fitness scores

cd ~/data/yeast/vcf/group

for group in $(ls *.bcf | sed 's/^1011Matrix\.//' | sed 's/\.bcf$//')
do
    echo "==>${group}"

    bcftools view 1011Matrix.${group}.bcf -Ov \
        -R ../region/region.bed -o ${group}.vcf
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
        > ../region/${group}.tsv
done
```

## Statistical analysis

### Combine vcf with fitness

- Add fitness data to vcf

```bash
mkdir ~/data/yeast/vcf/fit
mkdir ~/data/yeast/fitness
cd ~/data/yeast

for group in $(cat isolates/group.lst)
do
    echo "==> ${group}"

    tsv-join vcf/region/all.mut.vcf -k 1,2,3,4 \
        --filter-file <(cat vcf/region/${group}.tsv | tsv-filter --ne 5:0) \
        --append-fields 5,6,7 \
        > vcf/fit/${group}.fit.tsv
done

wc -l vcf/fit/*.fit.tsv
#  179 total
# 8341 muts, only 513 muts were existed among all wild groups

if [[ -f fitness/group.fit.tsv ]]; then
    rm fitness/group.fit.tsv
    echo 'Clear'
else
    echo 'Empty'
fi

# mv all fitness into one tsv file
for group in $(cat isolates/group.lst)
do
    echo "==> ${group}"
    
    cat vcf/fit/${group}.fit.tsv |
        awk -v col="${group}" '{print (col "\t" $0)}' \
        >> fitness/group.fit.tsv
done

cat fitness/group.fit.tsv | wc -l
#179

cat fitness/group.fit.tsv |
    tsv-uniq -f 2,3,4,5 |
    wc -l
#66
# totally 66 snps ( 1 snp may exists more than 1 group)
# but 1 snp may have different freq among groups

# other muts without occurring among all groups
cat vcf/region/all.mut.vcf |
    tsv-join -k 1,2,3,4 \
    --filter-file <(cat fitness/group.fit.tsv | tsv-select -f 2,3,4,5) \
    --exclude |
    awk '{print ("other" "\t" $0 "\t" "0" "\t" "-" "\t" "-")}' \
    > fitness/other.fit.tsv

# combine 2 parts
cat fitness/group.fit.tsv fitness/other.fit.tsv \
    > fitness/all.fit.tsv

cat fitness/all.fit.tsv | wc -l
#8454
# 8454 - (179 - 66) = 8341
# so some snps exist among groups

rm fitness/group.fit.tsv fitness/other.fit.tsv

# all nonsense mutations in wild
cat all.fit.tsv |
    tsv-filter --str-eq 7:Nonsense_mutation --str-ne 1:other \
    > nonsense.tsv
```

- Basic info

```bash
mkdir ~/data/yeast/fitness/freq
cd ~/data/yeast/fitness

# all repeated snps
cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-summarize -g 2,3,4,5 --count |
    tsv-filter --gt 5:1 |
    tsv-sort -r -nk 5,5 |
    parallel --col-sep "\t" -j 6 -k '
        gene=$(cat ../vcf/region/region_name.tsv |
                  tsv-filter --str-eq 3:{1} --le 4:{2} --ge 5:{2} |
                  tsv-select -f 1)
        echo -e "${gene}\t{1}\t{2}\t{3}\t{4}\t{5}"
    ' | mlr --itsv --omd cat

# number of snps among groups
cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-summarize -g 1 --count |
    tsv-sort -nk 2,2 -r |
    sed '1igroup\tnum' |
    mlr --itsv --omd cat

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

number of snps among groups:

| group        | num |
|--------------|-----|
| Wine         | 16  |
| Nature       | 14  |
| Beer         | 13  |
| Tree         | 12  |
| Distillery   | 11  |
| Unknown      | 10  |
| Fruit        | 10  |
| Clinical     | 10  |
| Soil         | 9   |
| Palm_wine    | 9   |
| Fermentation | 9   |
| Human        | 7   |
| Sake         | 6   |
| Industrial   | 6   |
| Bioethanol   | 6   |
| Insect       | 5   |
| Dairy        | 5   |
| Bakery       | 5   |
| Water        | 4   |
| Lab_strain   | 4   |
| Flower       | 4   |
| Cider        | 4   |

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

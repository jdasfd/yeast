# Check wild yeast vcf

The [article](https://www.nature.com/articles/s41586-022-04823-w) contains thousands of mutations from 21 genes (synonymous or non-synonymous mutations). The goal is to validate all mutations among wild yeasts.

`mkdir -p ~/data/yeast/info`

Download all supplementary files and tables from the [website](https://www.nature.com/articles/s41586-022-04823-w) and save them into info dir.

Steps followed here are processes to convert the mutation info into vcf format for better analyze.

There are some data from the [README.md](https://github.com/wang-q/pars#readme). Better completing those steps.

## Data preparation

### Split strains from 1002 genomes project according to ecological groups

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

### Split according to groups from the `1011Matrix.gvcf.gz`

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

### Extract all mutations in vcf-like format

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

### Locating mutations on the genome scale.

```bash
cd ~/data/yeast

# average fitness of muts were not all available
# remove those #DIV/0!
perl scripts/xlsx2csv.pl -f info/41586_2022_4823_MOESM9_ESM.xlsx \
    --sheet "Fig. 2abc" |
    sed '1d' |
    tsv-select -d, -f 1,3,2 |
    mlr --icsv --otsv cat |
    tr '-' '\t' |
    tsv-filter --not-iregex 5:# \
    > info/fit.tsv

for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"

    if ! [ -d gene/${gene} ]
    then
        mkdir gene/${gene}
    fi

    cat info/fit.tsv |
        grep "^$gene" \
        > gene/${gene}/${gene}.tsv
done

# extract 150 coding regions (detected)
for gene in $(cat gene/stdname.lst)
do
    echo "==> ${gene}"

    # count num of muts for every gene
    cat gene/${gene}/${gene}.tsv |
        tsv-select -f 1,2,3 |
        tsv-uniq |
        sort -nk 2,2 |
        wc -l

    echo ">${gene}" \
        > gene/${gene}/${gene}.mut.fa

    cat gene/${gene}/${gene}.tsv |
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
        -t gene/${gene}/${gene}.tsv \
        -a gene/${gene}/${gene}.aln.fa \
        > gene/${gene}/${gene}.mut.vcf
done
```

### Filtering with bcf of different groups using region

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

### Join analysis

- Join fitness data to vcf

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

rm fitness/group.fit.tsv

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
    tsv-summarize -g 7 --mean 6 --median 6 |
    mlr --itsv --omd cat

# other muts without occurring among all groups
cat vcf/region/all.mut.vcf |
    tsv-join -k 1,2,3,4 \
    --filter-file <(cat fitness/group.fit.tsv | tsv-select -f 2,3,4,5) \
    --exclude |
    awk '{print ("other" "\t" $0 "\t" "0" "\t" "-" "\t" "-")}' \
    > fitness/other.fit.tsv

cat fitness/other.fit.tsv |
    tsv-summarize -g 7 --mean 6 --median 6 |
    mlr --itsv --omd cat

# combine 2 parts
cat fitness/group.fit.tsv fitness/other.fit.tsv \
    > fitness/all.fit.tsv
```

`group.fit.tsv`:

| Nonsynonymous_mutation | 0.991225602326 | 0.99691229925  |
|------------------------|----------------|----------------|
| Synonymous_mutation    | 0.988202180851 | 0.991590123125 |
| Nonsense_mutation      | 0.944125370952 | 0.949334339282 |

`other.fit.tsv`:

| Nonsynonymous_mutation | 0.984944235527 | 0.98805154375  |
|------------------------|----------------|----------------|
| Synonymous_mutation    | 0.987772083901 | 0.988734612306 |
| Nonsense_mutation      | 0.934085408798 | 0.939873592125 |

- Count different fitness format and plot

```bash
cd ~/data/yeast/fitness

cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-summarize --group-by 7 --count |
    mlr --itsv --omd cat

cat all.fit.tsv |
    tsv-filter --str-ne 1:other |
    tsv-summarize --group-by 1,7 --mean 6 --median 6 |
    tsv-filter --str-ne 2:Nonsense_mutation |
    sed '1igroup\ttype\tmean\tmedian' \
    > fitness.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1])
p1 <- ggplot(fit, aes(x = group, y = mean, group = type)) +
geom_boxplot()
p2 <- ggplot(fit, aes(x = group, y = median, group = type)) +
geom_boxplot()
ggsave(p1, file = "mean.pdf")
ggsave(p2, file = "median.pdf")
' fitness.tsv
```

| Nonsynonymous_mutation | 89  |
|------------------------|-----|
| Synonymous_mutation    | 86  |
| Nonsense_mutation      | 4   |

| group        | type                   | mean           | median         |
|--------------|------------------------|----------------|----------------|
| Bakery       | Nonsynonymous_mutation | 0.993952946833 | 0.99691229925  |
| Bakery       | Synonymous_mutation    | 0.985707319827 | 0.985707319827 |
| Beer         | Nonsynonymous_mutation | 0.990627553393 | 0.9934689155   |
| Beer         | Synonymous_mutation    | 0.988379347109 | 0.989974881    |
| Bioethanol   | Nonsynonymous_mutation | 0.996937521125 | 0.996937521125 |
| Bioethanol   | Synonymous_mutation    | 0.990054031913 | 0.994400744    |
| Cider        | Nonsynonymous_mutation | 0.996937521125 | 0.996937521125 |
| Cider        | Synonymous_mutation    | 0.985707319827 | 0.985707319827 |
| Clinical     | Nonsynonymous_mutation | 0.993034458    | 0.995005548625 |
| Clinical     | Synonymous_mutation    | 0.990552478788 | 0.993249265125 |
| Dairy        | Nonsynonymous_mutation | 0.990048784583 | 0.9905381615   |
| Dairy        | Synonymous_mutation    | 0.995520676    | 0.995520676    |
| Distillery   | Synonymous_mutation    | 0.991897391109 | 0.9935399995   |
| Distillery   | Nonsynonymous_mutation | 0.99258668905  | 0.99691229925  |
| Fermentation | Synonymous_mutation    | 0.987350236817 | 0.99448643775  |
| Fermentation | Nonsynonymous_mutation | 0.990413197    | 0.99691229925  |
| Flower       | Nonsynonymous_mutation | 0.996937521125 | 0.996937521125 |
| Flower       | Synonymous_mutation    | 0.985707319827 | 0.985707319827 |
| Fruit        | Nonsynonymous_mutation | 0.993961673063 | 0.995394147    |
| Fruit        | Synonymous_mutation    | 0.984812363317 | 0.984694513    |
| Human        | Nonsynonymous_mutation | 0.996937521125 | 0.996937521125 |
| Human        | Synonymous_mutation    | 0.987182293481 | 0.9873332615   |
| Industrial   | Nonsynonymous_mutation | 0.993952946833 | 0.99691229925  |
| Industrial   | Synonymous_mutation    | 0.986496106134 | 0.98807367875  |
| Insect       | Nonsynonymous_mutation | 0.9882705025   | 0.99691229925  |
| Insect       | Synonymous_mutation    | 0.985707319827 | 0.985707319827 |
| Lab_strain   | Nonsynonymous_mutation | 0.996937521125 | 0.996937521125 |
| Lab_strain   | Synonymous_mutation    | 0.985707319827 | 0.985707319827 |
| Nature       | Nonsynonymous_mutation | 0.988430281375 | 0.989328105125 |
| Nature       | Synonymous_mutation    | 0.987996963359 | 0.989974881    |
| Palm_wine    | Nonsynonymous_mutation | 0.991533430125 | 0.99181591175  |
| Palm_wine    | Synonymous_mutation    | 0.989110202081 | 0.99119987325  |
| Sake         | Nonsynonymous_mutation | 0.993432876938 | 0.995394147    |
| Sake         | Synonymous_mutation    | 0.985707319827 | 0.985707319827 |
| Soil         | Nonsynonymous_mutation | 0.98570085225  | 0.99691229925  |
| Soil         | Synonymous_mutation    | 0.988496564976 | 0.991285810125 |
| Soil         | Nonsense_mutation      | 0.894531893404 | 0.894531893404 |
| Tree         | Nonsynonymous_mutation | 0.985150976651 | 0.984364901827 |
| Tree         | Synonymous_mutation    | 0.990710495678 | 0.994297732583 |
| Unknown      | Nonsynonymous_mutation | 0.9978449385   | 0.996937521125 |
| Unknown      | Nonsense_mutation      | 0.9937188485   | 0.9937188485   |
| Unknown      | Synonymous_mutation    | 0.989668183781 | 0.991980373    |
| Water        | Nonsynonymous_mutation | 0.996937521125 | 0.996937521125 |
| Water        | Synonymous_mutation    | 0.985707319827 | 0.985707319827 |
| Wine         | Nonsynonymous_mutation | 0.986167073736 | 0.990672412    |
| Wine         | Nonsense_mutation      | 0.9937188485   | 0.9937188485   |
| Wine         | Synonymous_mutation    | 0.985464730372 | 0.985222140917 |

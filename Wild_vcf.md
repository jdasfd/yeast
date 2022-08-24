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
```
# number of snps among groups
cat fitness/group.fit.tsv |
    tsv-summarize -g 1 --count |
    mlr --itsv --omd cat
```

`group.fit.tsv`:

| Nonsynonymous_mutation | 89  | 0.991225602326 | 0.99691229925  |
|------------------------|-----|----------------|----------------|
| Synonymous_mutation    | 86  | 0.988202180851 | 0.991590123125 |
| Nonsense_mutation      | 4   | 0.944125370952 | 0.949334339282 |

uniq all snps

| Nonsynonymous_mutation | 33  | 0.987204450667 | 0.98809565525  |
|------------------------|-----|----------------|----------------|
| Synonymous_mutation    | 30  | 0.990295101505 | 0.99229843675  |
| Nonsense_mutation      | 3   | 0.92759421177  | 0.904949830064 |

`other.fit.tsv`:

| Nonsynonymous_mutation | 6273 | 0.984944235527 | 0.98805154375  |
|------------------------|------|----------------|----------------|
| Synonymous_mutation    | 1836 | 0.987772083901 | 0.988734612306 |
| Nonsense_mutation      | 166  | 0.934085408798 | 0.939873592125 |

| Bakery       | 5   |
|--------------|-----|
| Beer         | 13  |
| Bioethanol   | 6   |
| Cider        | 4   |
| Clinical     | 10  |
| Dairy        | 5   |
| Distillery   | 11  |
| Fermentation | 9   |
| Flower       | 4   |
| Fruit        | 10  |
| Human        | 7   |
| Industrial   | 6   |
| Insect       | 5   |
| Lab_strain   | 4   |
| Nature       | 14  |
| Palm_wine    | 9   |
| Sake         | 6   |
| Soil         | 9   |
| Tree         | 12  |
| Unknown      | 10  |
| Water        | 4   |
| Wine         | 16  |

```bash
cd ~/data/yeast/fitness

# all nonsense mutations in wild
cat all.fit.tsv |
    tsv-filter --str-eq 7:Nonsense_mutation --str-ne 1:other \
    > nonsense.tsv

# 0.05 as cutoff
cat all.fit.tsv |
    tsv-filter --str-ge 8:0.05 --str-ne 7:Nonsense_mutation |
    tsv-select -f 1,6,7 |
    sed '1igroup\tfit\ttype' \
    > fit_ge.tsv

cat all.fit.tsv |
    tsv-filter --str-lt 8:0.05 --str-ne 7:Nonsense_mutation |
    tsv-select -f 1,6,7 |
    sed '1igroup\tfit\ttype' \
    > fit_lt.tsv

wc -l fit_*.tsv

#   106 fit_ge.tsv
#  8350 fit_lt.tsv
#  8456 total

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1])
p <- ggplot(fit, aes(x = group, y = fit, group = type, color = type)) +
      geom_point() +
      theme(axis.text.x = element_text(angle = 315))
ggsave(p, height = 6, width = 12, file = "fit_ge.pdf")
' fit_ge.tsv

Rscript -e '
library(ggplot2)
library(readr)
args <- commandArgs(T)
fit <- read_tsv(args[1])
p <- ggplot(fit, aes(x = group, y = fit, group = type, color = type)) +
      geom_point() +
      theme(axis.text.x = element_text(angle = 315))
ggsave(p, height = 6, width = 12, file = "fit_lt.pdf")
' fit_lt.tsv
```

| group        | type                   | num | mean           | median         |
|--------------|------------------------|-----|----------------|----------------|
| Bakery       | Nonsynonymous_mutation | 3   | 0.993952946833 | 0.99691229925  |
| Bakery       | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Beer         | Nonsynonymous_mutation | 4   | 0.996479855938 | 0.996937521125 |
| Beer         | Synonymous_mutation    | 4   | 0.987841100413 | 0.989974881    |
| Bioethanol   | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Bioethanol   | Synonymous_mutation    | 1   | 0.9984248515   | 0.9984248515   |
| Cider        | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Cider        | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Clinical     | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Clinical     | Synonymous_mutation    | 3   | 0.986496106134 | 0.98807367875  |
| Dairy        | Nonsynonymous_mutation | 2   | 0.989804096125 | 0.989804096125 |
| Dairy        | Synonymous_mutation    | 2   | 0.995520676    | 0.995520676    |
| Distillery   | Synonymous_mutation    | 3   | 0.989098255134 | 0.99588012575  |
| Distillery   | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Fermentation | Synonymous_mutation    | 4   | 0.991674809788 | 0.9964791265   |
| Fermentation | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Flower       | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Flower       | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Fruit        | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Fruit        | Synonymous_mutation    | 3   | 0.984050630801 | 0.98073725275  |
| Human        | Synonymous_mutation    | 4   | 0.984871288476 | 0.984035257125 |
| Human        | Nonsynonymous_mutation | 1   | 0.99691229925  | 0.99691229925  |
| Industrial   | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Industrial   | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Insect       | Nonsynonymous_mutation | 3   | 0.9882705025   | 0.99691229925  |
| Insect       | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Lab_strain   | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Lab_strain   | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Nature       | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Nature       | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Palm_wine    | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Palm_wine    | Synonymous_mutation    | 3   | 0.987538170968 | 0.99119987325  |
| Sake         | Nonsynonymous_mutation | 3   | 0.995917012333 | 0.99691229925  |
| Sake         | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Soil         | Synonymous_mutation    | 4   | 0.988496564976 | 0.991285810125 |
| Soil         | Nonsense_mutation      | 2   | 0.894531893404 | 0.894531893404 |
| Soil         | Nonsynonymous_mutation | 2   | 0.980069906875 | 0.980069906875 |
| Tree         | Nonsynonymous_mutation | 3   | 0.9882705025   | 0.99691229925  |
| Tree         | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Unknown      | Nonsynonymous_mutation | 3   | 0.99821723375  | 0.996962743    |
| Unknown      | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Water        | Nonsynonymous_mutation | 2   | 0.996937521125 | 0.996937521125 |
| Water        | Synonymous_mutation    | 2   | 0.985707319827 | 0.985707319827 |
| Wine         | Nonsynonymous_mutation | 1   | 0.99691229925  | 0.99691229925  |
| Wine         | Synonymous_mutation    | 1   | 0.9984248515   | 0.9984248515   |

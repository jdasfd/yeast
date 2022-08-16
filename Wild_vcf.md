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
```

- Count different fitness format and plot

```bash
cd ~/data/yeast/fitness


```

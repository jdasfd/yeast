#!/usr/bin/env bash

usage () { echo "bash vcf_num.sh <mut.snp.tsv>" 1>&2; exit; }
[ $# -lt 1 ] && usage

TSV_FILE=$1

if [ ! -e ${TSV_FILE} ]; then
    echo 1>&2 "[${TSV_FILE}] is not a tsv file" 1>&2;
    echo 1>&2 "Please provide CHROM, POS, REF, ALT, AF, AC, AN in tsv format" 1>&2;
    exit;
fi

echo "==> All SNPs"

all1=$(cat ${TSV_FILE} |
          tsv-summarize -g 1,2 --count |
          tsv-summarize -g 3 --count |
          tr '\t' '*' |
          paste -sd+ |
          bc)

all2=$(cat ${TSV_FILE} | wc -l)

if [[ $all1 != $all2 ]]; then
    echo 1&2 "Something wrong, please check manually!" 1>&2; exit;
else
    echo "Check OK!"
    echo "Total SNP number is ${all1}"
fi

echo
echo "==> Process all pos with SNP"

pos=$(cat ${TSV_FILE} |
          tsv-summarize -g 1,2 --count |
          tsv-summarize -g 3 --count |
          tsv-select -f 2 |
          paste -sd+ |
          bc)

echo "The number of SNP sites is ${pos}"
echo

echo "==> Process pos with 1 SNP"

one=$(cat ${TSV_FILE} |
          tsv-summarize -g 1,2 --count |
          tsv-summarize -g 3 --count |
          grep '^1' |
          tsv-select -f 2)

echo "The number of sites with 1 SNP is ${one}"
echo

echo "==> Process pos with 2 SNPs"

two=$(cat ${TSV_FILE} |
          tsv-summarize -g 1,2 --count |
          tsv-summarize -g 3 --count |
          grep '^2' |
          tsv-select -f 2)

echo "The number of sites with 2 SNPs is ${two}"
echo

echo "==> Process pos with 3 SNPs"

three=$(cat ${TSV_FILE} |
            tsv-summarize -g 1,2 --count |
            tsv-summarize -g 3 --count |
            grep '^3' |
            tsv-select -f 2)

echo "The number of sites with 3 SNPs is ${three}"
echo

echo "==> Final result"

echo -e "| Mut | Num |"
echo -e "| --- | --- |"
echo -e "| All | ${all1} |"
echo -e "| Pos | ${pos} |"
echo -e "| One_SNP | ${one} |"
echo -e "| Two_SNPs | ${two} |"
echo -e "| Three_SNPs | ${three} |"

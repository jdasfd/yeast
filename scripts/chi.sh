cd ~/data/yeast

# genome length
genom_len=$(cat GENOMES/S288c/chr.sizes |
                tsv-filter --str-ne 1:Mito |
                tsv-select -f 2 |
                paste -sd+ |
                bc)

# gene length
gene_len=$(spanr stat GENOMES/S288c/chr.sizes ../mrna-structure/gene-filter/genes.merge.yml |
              tsv-filter -H -d, --str-ne chr:all |
              tsv-select -d, -f 4 |
              sed 1d |
              paste -sd+ |
              bc)

# other region
other_len=$(echo ${genom_len}-${gene_len} | bc)

# 21 gene length
gene_21_len=$(faops size gene/gene.mut.fa |
                  tsv-select -f 2 |
                  paste -sd+ |
                  bc)

# genome snps
genom_vcf=$(cat vcf/all.snp.tsv | tsv-uniq -f 1,2 | wc -l)
gene_vcf=$(cat vcf/gene_all.snp.tsv | tsv-uniq -f 1,2 | wc -l)
not_vcf=$(cat vcf/not_gene.snp.tsv | tsv-uniq -f 1,2 | wc -l)
gene_21_vcf=$(cat vcf/gene_21.snp.tsv | tsv-uniq -f 1,2 | wc -l)

echo "==> genome vs gene"
echo -e "$genom_len\t$gene_len\t$genom_vcf\t$gene_vcf" |
    parallel --colsep '\t' -j 1 -k '
        Rscript -e "
            x <- matrix(c({1},{3},{2},{4}), ncol = 2)
            old.warn <- options()$warn
            options(warn = -1)
            x
            chisq.test(x)
        "
'

echo "==> gene vs other"
echo -e "$gene_len\t$other_len\t$gene_vcf\t$not_vcf" |
    parallel --colsep '\t' -j 1 -k '
        Rscript -e "
            x <- matrix(c({1},{3},{2},{4}), ncol = 2)
            old.warn <- options()$warn
            options(warn = -1)
            x
            chisq.test(x)
        "
'

echo "==> gene vs 21"
echo -e "$gene_len\t$gene_21_len\t$gene_vcf\t$gene_21_vcf" |
    parallel --colsep '\t' -j 1 -k '
        Rscript -e "
            x <- matrix(c({1},{3},{2},{4}), ncol = 2)
            old.warn <- options()$warn
            options(warn = -1)
            x
            chisq.test(x)
        "
'

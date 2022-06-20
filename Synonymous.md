# Synonymous mutation among yeast

This markdown records the process of repeating the analysis from this [article](https://www.nature.com/articles/s41586-022-04823-w).

## Seq files

```bash
mkdir -p /mnt/e/data/yeast/ena
cd /mnt/e/data/yeast/ena

# download SraRunTable.txt to this doc
# PRJNA750109, search it in NCBI SRA Run Selector

# SraRunTable.txt is in .csv format
cat SraRunTable.txt |
    mlr --icsv --otsv cat |
    tsv-select -H -f Experiment,"Sample\ Name",Bases \
    > SraRunTable.tsv

cat SraRunTable.tsv |
    sed '1 s/^/#/' |
    keep-header -- tsv-sort -k2,2 -k3,3nr |
    tsv-uniq -H -f "Sample\ Name" --max 1 |
    mlr --itsv --ocsv cat \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv
# viewing md format of each sequencing file

aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt
# aria2c could be used, although it was better using aspera (shell need to be modified if using aspera)

md5sum --check ena_info.md5.txt
# check if there was any mistake during downloading
```
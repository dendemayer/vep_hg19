#!/usr/bin/env bash

# this scripts takes 
file=$1
wget --no-clobber --directory-prefix=ref_seq https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/$file
file=ref_seq/$file
echo file to download is $file
MD5SUM=$(md5sum $file | cut -d' ' -f1)
file_in_MD5SUMS=$(grep $MD5SUM ref_seq/MD5SUMS | cut -d' ' -f3)
#echo file_in_MD5SUMS is $file_in_MD5SUMS
if [[ '$file'=='$file_in_MD5SUMS' ]]; then
    echo succesfully downloaded file $file with md5sum $MD5SUM
else
    echo file: $file and file_in_MD5SUMS: $file_in_MD5SUMS .
    echo something went wrong while download of $file
    exit 1
fi

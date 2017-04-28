#!/bin/bash
schema=$1;
ncbi=$2;
source=$3;
out=$schema'@'$ncbi'@'$source;
/usr/bin/mysql -h 172.17.100.150 -pserver <query.sql |sed -e 's/[[:space:]]*$//g' |sed '1d' | sed -e 's/\\n/\n/' -e 's/\*$//' >$out.fa

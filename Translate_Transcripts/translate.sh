#!/bin/bash
cd `pwd`
#USAGE:  $1= Speciesname $2= Trintyassembly  $3= blastdb $4= threads $5= amino acid for trans

mkdir $1_Translation
cd $1_Translation
time blastpath/tblastx -num_threads $4 -query $2 -db $3 -outfmt 6 -evalue 1e-10 -out $1.tblastx
while [ ! -f "$1.tblastx" ];
do
	sleep 15
done

export WISECONFIGDIR=basedir/wise2.2.0/wisecfg/
perl basedir/TranslationsPipeline.pl $1.tblastx $2 $1 $5 $4

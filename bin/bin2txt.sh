#!/bin/sh
#
echo $1
echo $2
echo $3
echo $4
echo $5
echo "for PNMHD, like: sh bin2txt.sh nx ny dump f 9" 
#echo "for GRMHD, like: sh bin2txt.sh nx ny dump d 26" 
echo "for GRMHD, like: sh bin2txt.sh nx ny dump d 34" 
for fil in `/bin/ls $3????.dat.bin`
do 
	prefa=`echo $fil | sed "s/\./ /"`
	pref=`echo $prefa | awk '{print $1}'`
	echo $pref
	bin2txt 1 2 0 0 3 $1 $2 1 1 $pref.dat.bin $pref.dat.txt $4 $5
	cat $pref.dat.bin.head $pref.dat.txt > $pref.dat
        rm $pref.dat.txt
	echo $pref done
done

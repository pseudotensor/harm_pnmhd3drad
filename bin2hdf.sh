#!/bin/sh
#
for fil in `/bin/ls dump????.dat`
do 
	prefa=`echo $fil | sed "s/\./ /"`
	pref=`echo $prefa | awk '{print $1}'`
	echo $pref
	./bin2txt 1 3 0 128 128 128 $pref.dat $pref.hdf f 7 
	echo $pref done
done

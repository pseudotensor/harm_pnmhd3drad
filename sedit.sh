#!/bin/bash

#sed -e "s/HOR=$3/HOR=$6/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1

# start out at 32,32, specific HOR
ovar1=32
ovar0=32

var1=32
while [ "$var1" != 512 ]
do
var0=32
while [ "$var0" != 512 ]
do
  echo "$ovar0 $ovar1 -> $var0 $var1 "
#sed -e "s/#define N1 $ovar0/#define N1 $var0/" global.h > global.h.1 ; cp global.h.1 global.h ; rm global.h.1
#sed -e "s/#define N2 $ovar1/#define N2 $var1/" global.h > global.h.1 ; cp global.h.1 global.h ; rm global.h.1
	make
  mv ./bin/mhd3d ./bin/mhd3d.$var0-$var1


  ovar0=$((var0))
  ovar1=$((var1))
  var0=$(($var0*2))
done
  var1=$(($var1*2))
done

exit 0


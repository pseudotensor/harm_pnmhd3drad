#!/bin/bash
# start with H/R=.1
  make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd1281
  qsub ./q/thindisk.pt

sed -e "s/HOR=0.1/HOR=0.05/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
sed -e "s/acoef=56.1937/acoef=137.767/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
sed -e "s/td1281/td12805/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
  make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd12805
  qsub ./q/thindisk.pt

sed -e "s/HOR=0.05/HOR=0.025/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
sed -e "s/acoef=137.767/acoef=315.941/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
sed -e "s/td12805/td128025/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
  make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd128025
  qsub ./q/thindisk.pt

sed -e "s/HOR=0.025/HOR=0.0125/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
sed -e "s/acoef=315.941/acoef=707.359/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
sed -e "s/td128025/td1280125/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
	make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd1280125
  qsub ./q/thindisk.pt

# these last 2 are too expensive
#sed -e "s/HOR=0.0125/HOR=0.00625/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
#sed -e "s/acoef=707.359/acoef=1551.66/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
#sed -e "s/td1280125/td12800625/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
#	make -f makefile.pt
#  mv ./bin/mhd3d ./bin/mhd3dtd12800625
#  qsub ./q/thindisk.pt
#
#sed -e "s/HOR=0.00625/HOR=0.003125/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
#sed -e "s/acoef=1551.66/acoef=3355.35/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
#sed -e "s/td12800625/td128003125/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
#	make -f makefile.pt
#  mv ./bin/mhd3d ./bin/mhd3dtd128003125
#  qsub ./q/thindisk.pt

exit 0


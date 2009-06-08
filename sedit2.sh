#!/bin/bash

#  make -f makefile.pt
#  mv ./bin/mhd3d ./bin/mhd3dtd1
#  qsub ./q/thindisk.pt

#sed -e "s/HOR=0.1/HOR=0.05/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
#sed -e "s/acoef=56.1937/acoef=137.767/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
#sed -e "s/td1/td05/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
#	make -f makefile.pt
#  mv ./bin/mhd3d ./bin/mhd3dtd05
#  qsub ./q/thindisk.pt

sed -e "s/HOR=0.05/HOR=0.025/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
sed -e "s/acoef=137.767/acoef=315.941/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
sed -e "s/td05/td025/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
  make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd025
  qsub ./q/thindisk.pt

sed -e "s/HOR=0.025/HOR=0.0125/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
sed -e "s/acoef=315.941/acoef=707.359/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
sed -e "s/td025/td0125/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
	make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd0125
  qsub ./q/thindisk.pt

sed -e "s/HOR=0.0125/HOR=0.00625/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
sed -e "s/acoef=707.359/acoef=1551.66/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
sed -e "s/td0125/td00625/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
	make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd00625
  qsub ./q/thindisk.pt

sed -e "s/HOR=0.00625/HOR=0.003125/" analsol.c > analsol.c.1 ; cp analsol.c.1 analsol.c ; rm analsol.c.1
sed -e "s/acoef=1551.66/acoef=3355.35/" init.c > init.c.1 ; cp init.c.1 init.c ; rm init.c.1
sed -e "s/td00625/td003125/" ./q/thindisk.pt > ./q/thindisk.pt.1 ; cp ./q/thindisk.pt.1 ./q/thindisk.pt ; rm ./q/thindisk.pt.1
	make -f makefile.pt
  mv ./bin/mhd3d ./bin/mhd3dtd003125
  qsub ./q/thindisk.pt


exit 0


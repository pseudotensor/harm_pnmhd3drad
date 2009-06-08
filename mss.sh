#!/bin/sh
#
for fil in `/bin/ls`
do
    kftp -i -n <<HERE
         open mss.ncsa.uiuc.edu
         user jmckinne
         put $fil
         chmod 0644 $fil
         close
         quit
HERE
done


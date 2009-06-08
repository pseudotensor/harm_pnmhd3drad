#!/bin/sh
sed 's/POSTPROC 1/POSTPROC 0/' global.h > globaltemp.h
cp globaltemp.h global.h

#!/bin/sh
sed 's/POSTPROC 0/POSTPROC 1/' global.h > globaltemp.h
cp globaltemp.h global.h

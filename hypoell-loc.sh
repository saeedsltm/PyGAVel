#!/bin/bash
# SCRIPT FOR RUNNING HYPOELLIPSE

# PHASE FILE 
filepick=hypoel.pha

# VELOCITY FILE
filevel=${1}.prm

# STATION FILE
filesta=hypoel.sta

# $root.out : FINAL LOCATION FILE
# $root.sum : SUMMARY FILE
# $root.arc : POLARITIES
root=${1}

# PARAMETER FILE
filepar=default.cfg

# REFERNCE TIME
date=19800101

# CONTROL FILE
filecom=${1}.ctl
echo "stdin">$filecom
echo "y">>$filecom
echo $root>>$filecom
echo "stdout">>$filecom
echo "screen">>$filecom
echo "y">>$filecom
echo $root.sum>>$filecom
echo "y">>$filecom
echo $root.arc>>$filecom
echo "n">>$filecom
echo "n">>$filecom

# JUMT TO PARAMETER FILE
echo "jump" $filepar>>$filecom

# JUMT TO VELOCITY FILE
echo "jump" $filevel>>$filecom

# USED STATIONS
echo "begin station list +1" $date >>$filecom

# JUMT TO STATION FILE
echo "jump" $filesta >>$filecom

# USED PHASE
echo "arrival times next">>$filecom

# JUMT TO PHASE FILE
echo "jump" $filepick >>$filecom

# RUN HYPOELLIPSE
echo "hymain < $filecom 2>&1 >/dev/null" > ${1}_runHE

chmod +x ${1}_runHE

./${1}_runHE > $root.out
rm ./${1}_runHE


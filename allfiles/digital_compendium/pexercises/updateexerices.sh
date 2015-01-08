#!/bin/bash

SVNDIR=~/Documents/tkt4140/Pythonovingar

#echo $SVNDIR/tkt4140/Pythonovingar/poving?

svndir_length=${#SVNDIR}

WD=$(pwd)

#echo 'you are here:' $WD


for oving in $(ls -d $SVNDIR/poving?) ; do
    local_oving=${oving:svndir_length+1}
    echo "mkdir" $local_oving
    mkdir $local_oving

    echo "cp" $oving/index.html  $local_oving
    cp $oving/index.html  $local_oving

    for files in $(ls $oving/*.png); do
	echo "cp" $files $local_oving
	cp $files $local_oving
    done

done


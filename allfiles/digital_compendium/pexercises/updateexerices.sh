<<<<<<< HEAD
#!/bin/bash

#SVNDIR=~/Documents/tkt4140/Pythonovingar
SVNDIR=~/Dokument/tkt4140/Pythonovingar

#echo $SVNDIR/tkt4140/Pythonovingar/poving?

svndir_length=${#SVNDIR}

WD=$(pwd)

echo 'you are here:' $WD


for exercise in $(ls -d $SVNDIR/pexercise?) ; do
    local_exercise=${exercise:svndir_length+1}
    echo "mkdir" $local_exercise
    mkdir $local_exercise

    echo "cp" $exercise/index.html  $local_exercise
    cp $exercise/index.html  $local_exercise

    for files in $(ls $exercise/*.png); do
	echo "cp" $files $local_exercise
	cp $files $local_exercise
    done

done

=======
#!/bin/bash

#SVNDIR=~/Documents/tkt4140/Pythonovingar
SVNDIR=~/Dokument/tkt4140/Pythonovingar

#echo $SVNDIR/tkt4140/Pythonovingar/poving?

svndir_length=${#SVNDIR}

WD=$(pwd)

echo 'you are here:' $WD


for exercise in $(ls -d $SVNDIR/pexercise?) ; do
    local_exercise=${exercise:svndir_length+1}
    echo "mkdir" $local_exercise
    mkdir $local_exercise

    echo "cp" $exercise/index.html  $local_exercise
    cp $exercise/index.html  $local_exercise

    for files in $(ls $exercise/*.png); do
	echo "cp" $files $local_exercise
	cp $files $local_exercise
    done

done

>>>>>>> master

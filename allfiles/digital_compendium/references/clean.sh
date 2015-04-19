<<<<<<< HEAD
#!/bin/bash

rm papers.pub

for file in *.bib
do
    echo $file
    publish import $file

done

=======
#!/bin/bash

rm papers.pub

for file in *.bib
do
    echo $file
    publish import $file

done

publish export papers.bib
>>>>>>> master

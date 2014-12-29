#!/bin/bash

rm papers.pub

for file in *.bib
do
    echo $file
    publish import $file

done


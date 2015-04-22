for a in `ls *.pdf`; do
    b=${a%.pdf}
    convert -density 600 -depth 8 -quality 85 -transparent black $b.pdf $b.png
done

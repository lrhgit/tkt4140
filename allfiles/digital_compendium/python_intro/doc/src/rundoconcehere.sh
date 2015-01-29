#!/bin/sh
name=bumpy
wrap=$name


# Plain HTML
doconce format html $wrap --html_style=bootstrap_bluegray --html_links_in_new_window --device=screen 
if [ $? -ne 0 ]; then echo "doconce could not compile document"; exit; fi
doconce replace '<hr class="figure">' ' ' *.html
doconce replace 'https://raw.github.com/hplgit/doconce/master/bundled/html_styles/style_bootstrap/css/bootstrap_bluegray.css' 'bootstrap_bluegray.css' *.html
sed -i '1i <!DOCTYPE HTML>' *.html
doconce split_html $wrap.html --pagination

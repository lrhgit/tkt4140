#!/bin/sh
name=main
wrap=$name

doconce format html $wrap --html_style=bootstrap_bluegray --examples_as_exercises --html_links_in_new_window --device=screen 
doconce replace '<hr class="figure">' ' ' *.html
doconce replace 'https://raw.github.com/hplgit/doconce/master/bundled/html_styles/style_bootstrap/css/bootstrap_bluegray.css' 'bootstrap_bluegray.css' *.html
sed -i '1i <!DOCTYPE HTML>' *.html
doconce split_html $wrap.html --pagination

#--urlcheck
#--section_numbering=on (problem with example as exercises)
# --html_bootstrap_navbar=off turns off the navigation bar

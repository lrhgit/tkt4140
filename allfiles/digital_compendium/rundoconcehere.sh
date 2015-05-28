#!/bin/sh
set -x # display shell commands prior to execution
name=main
wrap=$name

doconce format html $wrap --html_style=bootstrap_bluegray --examples_as_exercises --html_links_in_new_window --device=screen  --html_figure_hrule=none --html_DOCTYPE --encoding=utf-8
doconce split_html $wrap.html --pagination --nav_button=top+bottom
#exit

# Wise to run latex now and then to get error messages in equations
# (will lead to display errors in MathJax anyway...)
doconce format pdflatex $wrap --device=paper --encoding=utf-8
doconce ptex2tex $wrap envir=Verbatim  # very simple code envir
pdflatex $wrap
bibtex $wrap
pdflatex $wrap
pdflatex $wrap

#--urlcheck
#--section_numbering=on (problem with example as exercises)
# --html_bootstrap_navbar=off turns off the navigation bar

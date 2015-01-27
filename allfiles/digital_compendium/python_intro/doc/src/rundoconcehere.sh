#!/bin/sh
name=bumpy
wrap=$name


# Plain HTML
#doconce format html $wrap --html_style=bootstrap_bluegray --html_links_in_new_window --device=screen 
doconce format html $wrap --html_style=bootstrap_bluegray --examples_as_exercises --html_links_in_new_window --device=screen  --html_figure_hrule=none --html_DOCTYPE --encoding=utf-8
if [ $? -ne 0 ]; then echo "doconce could not compile document"; exit; fi
doconce split_html $wrap.html --pagination --nav_button=top+bottom


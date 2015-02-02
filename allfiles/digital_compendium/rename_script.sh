#!/bin/sh

# This script replaces text string in all doconce source files and python programs in the project. Intended for changing strings when changing folder structure etc.

oldname='ode_schemes'
newname='ODEschemes'

# change in all doconce-files
find . -name '*.do.txt' -print -exec doconce replace $oldname $newname {} \;

# change in all python files
find . -name '*.py' -print -exec doconce replace $oldname $newname {} \;

# remove new ~~files
find . -name '*~~' -ctime -1 -print -exec rm {} \;
find . -name '*~' -ctime -1 -print -exec rm {} \;

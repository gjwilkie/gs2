#!/bin/bash
#A small script to extract the input file from the new netcdf diagnostics file.

#File to look in
FIL=${1:-"*.cdf"}

#Here we:
#1: Get the input_file variable from the netcdf file
#2: Only print lines between '${VAR} = "' and '" ;' (i.e. ignore header and footer)
#3: Convert \\n to new lines
#4: Delete empty lines
#5: Ignore first line
#6: Ignore last line  
#7: Fix " style quotes
#8: Fix ' style quotes
ncdump -v input_file ${FIL} |\
 sed -n '/input_file = /,/" ;/p' |\
 sed 's|\\\\n|\n|g' |\
 sed '/^ *$/d' |\
 tail -n+2 |\
 head -n-2 |\
 sed 's|\\\"|\"|g' |\
 sed "s|\\\'|\'|g"

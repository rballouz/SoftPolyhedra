#!/bin/bash

#-e flag enables use of backslashes
echo -e "#declare $1 = mesh2 {"
echo -e "\tvertex_vectors {"

#Get total number of lines from faces file
nrecs=`awk 'END {print NR}' $1`
awk -v nrecs=$nrecs 'BEGIN {print "\t"nrecs*3","}\
                    {if (NR<nrecs) print "\t<"$1","$2","$3"> , <"$4","$5","$6"> , <"$7","$8","$9">,";\
                    else print "\t<"$1","$2","$3"> , <"$4","$5","$6"> , <"$7","$8","$9"> "}' $1

echo -e "\t}"

echo -e "\tface_indices {"
awk -v nrecs=$nrecs 'BEGIN {print "\t"nrecs","; i=0}\
                    {if (NR<nrecs) print "\t<"i","i+1","i+2">,";\
                    else print "\t<"i","i+1","i+2">";\
                    i=i+3}' $1
echo -e "\t}"

echo -e "\tpigment {Yellow}"
echo -e "}"
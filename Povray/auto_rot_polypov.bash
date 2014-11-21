#!/bin/bash

rm *.png
i=0
while [ $i -lt 360 ]
do
    sed "s/z\*0/z\*$((i+1))/" poly_demo.pov > poly_demo_def.pov
    povray Output_File_Name=poly_demo_`printf %04d $((i+1))`.png poly_demo_def.pov
    i=$((i+1))
done


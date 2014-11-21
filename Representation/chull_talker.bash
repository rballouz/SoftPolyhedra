#!/bin/bash
#usage: ./chull_talker.bash

for i in `seq 7`
do
    j=$((i+5))
    echo $j
    python random_sphere.py $j > hull_3d.in
    ./chull_3d_flt > polyfaces_`printf %03d $j`.dat
done

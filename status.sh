#!/bin/bash
dir="./jobs"
for i in $(ls $dir/*.log | sort -t "-" -n -k4,4 -k3,3 -k2,2 -k1,1); do
    printf "%3d %s\n" "$(cat $i | grep ^mean | tail -n +1 | wc -l)" "$i"
done

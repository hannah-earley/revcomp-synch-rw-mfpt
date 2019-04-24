#!/bin/bash
dir="./jobs"
for i in $(ls $dir/*.log | sort -t "-" -n -k1,1 -k2,2 -k3,3 -k4,4); do
    printf "%3d %s\n" "$(cat $i | grep ^mean | tail -n +1 | wc -l)" "$i"
done

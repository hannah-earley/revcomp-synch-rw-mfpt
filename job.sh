#!/bin/bash
bias=1
dist=1
widt=1
n=1000
m=1000
s=1000

while getopts "b:d:w:n:m:s:x:" o; do
    case "${o}" in
        b)
            bias="${OPTARG}"
            ;;
        d)
            dist="${OPTARG}"
            ;;
        w)
            widt="${OPTARG}"
            ;;
        n)
            n="${OPTARG}"
            ;;
        m)
            m="${OPTARG}"
            ;;
        s)
            s="${OPTARG}"
            ;;
        x)
            m="${OPTARG}"
            s="${OPTARG}"
            ;;
    esac
done

dir="./jobs"
file="2d-$bias-$widt-$dist"
log="$dir/$file.log"
dat="$dir/$file.dat"

mkdir -p "$dir"
echo "#" ./walk -2 -v -b "$bias" -d "$dist" -w "$widt" \
     -n "$n" -m "$m" -s "$s" -p "$dat" >> "$log"
./walk -2 -v -b "$bias" -d "$dist" -w "$widt" -n "$n" \
       -m "$m" -s "$s" -p "$dat" 2>&1 | tee -a "$log"

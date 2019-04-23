#!/bin/bash
bias=1
dist=1
widt=1
n=1000
m=1000
s=1000

usage () {
    echo "Usage: $0"
    echo "         [-b bias] [-d distance] [-w width] [-n n_wlkrs]"
    echo "         [-m n_meas] [-s smpl_wndw] [-x n_meas^smpl_wndw]"
}

while getopts ":b:d:w:n:m:s:x:h" o; do
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
        *|h)
            usage
            exit 1
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

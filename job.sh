#!/bin/bash
bias=0
dist=1
widt=1
suff=

usage () {
    echo "Usage: $0 [OPTIONS]"
    echo "  ADDITIONAL OPTIONS:"
    echo "    -S suffix    suffix/variant for the filename/output,"
    echo "                 useful for concurrent runs of a given job"
    echo
    echo "  - Otherwise options same as for ./walk"
    echo "  - -2v is added to options automatically"
    echo "  - Output filenames will be generated from"
    echo "    the supplied bias and supplied to ./walk"
    echo
    ./walk -h
}

options () {
    while getopts ":b:d:w:S:h" o; do
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
            S)
                suff="-${OPTARG}"
                ;;
            h)
                usage
                exit 1
                ;;
            \?|:)
                ;;
        esac
    done
}

options "$@"
while [ "$OPTIND" -lt "$#" ]; do
    # keep processing arguments after unexpected options...
    OPTIND+=1
    options "$@"
done

dir="./jobs"
file="2d-$bias-$widt-$dist$suff"
log="$dir/$file.log"
csv="$dir/$file.csv"
dat="$dir/$file.dat"

mkdir -p "$dir"
echo "#" "$0" "$@" >> "$log"
echo "#" ./walk "$@" -2 -v -p "$dat" -q "$csv" >> "$log"
./walk "$@" -2 -v -p "$dat" -q "$csv" 2>&1 | tee -a "$log"

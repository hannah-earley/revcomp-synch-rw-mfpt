#!/bin/bash
bias=0
dist=1
widt=1
suff=
disp_path=0

usage () {
    echo "Usage: $0 [OPTIONS]"
    echo "  ADDITIONAL OPTIONS:"
    echo "    -S suffix    suffix/variant for the filename/output,"
    echo "                 useful for concurrent runs of a given job"
    echo "    -P           display generated path and exit"
    echo
    echo "  - Otherwise options same as for ./walk"
    echo "  - -2v is added to options automatically"
    echo "  - Output filenames will be generated from"
    echo "    the supplied bias and supplied to ./walk"
    echo
    ./walk -h
}

argv=()
options () {
    while getopts ":b:d:w:S:hP" o; do
        case "${o}" in
            b)
                bias="${OPTARG}"
                argv+=(-b "$bias")
                ;;
            d)
                dist="${OPTARG}"
                argv+=(-d "$dist")
                ;;
            w)
                widt="${OPTARG}"
                argv+=(-w "$widt")
                ;;
            S)
                suff="-${OPTARG}"
                ;;
            h)
                usage
                exit 1
                ;;
            \?|:)
                argv+=("-${OPTARG}")
                ;;
            P)
                disp_path=1
                ;;
        esac
    done
}

options "$@"
while [ "$OPTIND" -le "$#" ]; do
    # keep processing arguments after unexpected options...
    argv+=("${!OPTIND}")
    OPTIND+=1
    options "$@"
done

dir="./jobs"
file="2d-$bias-$widt-$dist$suff"
log="$dir/$file.log"
csv="$dir/$file.csv"
dat="$dir/$file.dat"

if [ "$disp_path" -eq "1" ]; then
    echo "$dir/$file"
    exit 0
fi

argv+=(-2 -v)
argv+=(-p "$dat")
argv+=(-q "$csv")

mkdir -p "$dir"
echo "#" "$0" "$@" >> "$log"
echo "#" ./walk "${argv[@]}" >> "$log"

set -o pipefail
./walk "${argv[@]}" 2>&1 | tee -a "$log"
#!/usr/bin/env bash

function usage () {
	echo ""
    echo "unzip *.tar.bz2 file and then gzip its contents"
    echo ""
    echo "usage: bash bunzip_gzip_file.sh <directory.tar.bz2>"
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

f = $1

tar -vxjf $f
file="${f%\.tar*}"
gzip $file/*

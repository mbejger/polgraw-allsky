#!/bin/bash

find . -name "$1*.sum" -exec cat {} >> $1.summary \;
awk '{if($6) print $0}' $1.summary|wc

exit 0

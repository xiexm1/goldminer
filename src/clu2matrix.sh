#!/usr/bin/bash

set -euox pipefail

inputfile=$1
prefix=$2
clupath=$3
output=`dirname $inputfile`

clu2tb $inputfile $prefix $output

tb2matrix $prefix $output

sed -i -e 's/clu://g' -e 's/|/,/g' $output/${prefix}.matrix

matrix_trans $prefix $output $clupath

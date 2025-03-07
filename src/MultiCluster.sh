#!/usr/bin/bash
info(){
    echo "Usage:[goldminer] MultiCluster [options] file"
    echo "Description:"
    echo "        Step4: Build clusters colinear network in all genomes (N ≥ 3)"
    echo "        Example:" 
    echo "        [goldminer] MultiCluster -d output_dir -i input_dir -t input_table.txt -p output_prefix"
    echo ""
    echo "Options:"
    echo "        -d DIR, --dir=DIR"
    echo "                directory of output files"
    echo ""
    echo "        -i INPUT, --input=INPUT"
    echo "                directory of input files"
    echo ""
    echo "        -c CLUSTER, --cluster=CLUSTER"
    echo "                directory of cluster files"
    echo ""
    echo "        -t TABLE, --table=TABLE"
    echo "                genome list file"
    echo ""
    echo "        -p PREFIX, --prefix=PREFIX"
    echo "                prefix of output file"
    echo ""
    echo "        -h, --help"
    echo "                print help message and exit"
    echo ""
}

export -f info

output_dir=""
input_dir=""
input_table=""
output_prefix=""

while getopts ":d:i:c:t:p:h" opt; do
    case $opt in
        d) output_dir="$OPTARG" ;;
        i) input_dir="$OPTARG" ;;
        c) clu_path="$OPTARG" ;;
        t) input_table="$OPTARG" ;;
        p) output_prefix="$OPTARG" ;;
        h) info; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done

shift $((OPTIND - 1))

# if [ $# -eq 0 ]; then
#     echo "Error: No file specified."
#     info
#     exit 1
# fi

if [ ! -d "$output_dir" ]; then
    echo "Error: Output directory $output_dir does not exist." >&2
    exit 1
fi
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory $input_dir does not exist." >&2
    exit 1
fi
if [ ! -f "$input_table" ]; then
    echo "Error: Input table file $input_table does not exist." >&2
    exit 1
fi

gmlst=$(tr '\n' ' ' < "$input_table")

for gmf in $gmlst; do
    for gmt in $gmlst; do
        file="$input_dir/$gmf.$gmt.link"
        if [ -f "$file" ]; then
            awk '$4 != 0 {print $1":"$2,$3":"$4,$5}' "$file" >> "$output_dir/$output_prefix.tmp"
        else
            echo "$gmf $gmt" >> "$output_dir/$output_prefix.wo.link"
        fi
    done
done

module load mcl/12-068 

mcl_i=2
ResultsDir="$output_dir/${output_prefix}_$(date +%Y%m%d)"
mkdir -p "$ResultsDir"

if [ -f "$output_dir/$output_prefix.tmp" ]; then
    mcl "$output_dir/$output_prefix.tmp" --abc -I "$mcl_i" -o "$ResultsDir/NetMatrix_${mcl_i}_${output_prefix}"
else
    echo "Error: Temporary file $output_dir/$output_prefix.tmp does not exist." >&2
    exit 1
fi

cd "$ResultsDir" || { echo "Error: Failed to change directory to $ResultsDir." >&2; exit 1; }

export PATH=$PATH:/data2/user2/xiexm/projs/GeneODL/GoldMiner/bin

clu2matrix "$ResultsDir/NetMatrix_${mcl_i}_${output_prefix}" "${mcl_i}_${output_prefix}" ${clu_path}
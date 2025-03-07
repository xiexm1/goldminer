#!/usr/bin/env python3
import sys

try:
    inputfile = sys.argv[1]
except IndexError:
    print("ERR:Invalid input file")
    sys.exit()

prefix = sys.argv[2]
outputfile = sys.argv[3]

tb_file = open(f"{outputfile}/{prefix}.csv", 'w')
with open(inputfile, "r") as f:
    num = 1
    for line in f:
        data_list = line.strip().split('\t')
        for data in data_list:
            subdata = data.split(':')
            try:
                int(subdata[1])
                subdata.append("HOC"+str(num).zfill(6))
                tb_file.write('\t'.join(subdata)+'\n')
            except ValueError:
                continue
        num += 1
tb_file.close()

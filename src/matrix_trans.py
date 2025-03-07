#!/usr/bin/env python3
# Author:Xie Xiaoming
# Func:For gene backbone matrix transform to size and gene matrix.
import sys
prefix = sys.argv[1]
output_path = sys.argv[2] + '/'
clu_path = sys.argv[3] + '/'

with open(output_path + prefix + ".matrix", 'r') as f:
    file = f.read()
    data = file.split("\n")

files = data[0].split()
clu_data = []
for file in files[1:]:
    with open(clu_path+file+'.clu') as f:
        clu_data.append(f.readlines())

size_csv_file = open(output_path + prefix + '.csv', 'w')
size_csv_file.write("genome\thoc\tgene_count\n")

head = "\t".join(data[0].strip().split()) + "\n"
size_file = open(output_path + prefix + '.size', 'w')
size_file.write(head)
clu_file = open(output_path + prefix + '.clu', 'w')
clu_file.write(head)
chr_file = open(output_path + prefix + '.chr', 'w')
chr_file.write(head)
chr_sum = open(output_path + prefix + '.chrsum', 'w')
chr_sum.write('Clu\tChr\n')

for subdata in data[1:]:
    backbone_node = subdata.strip().split()
    ind = 0
    s = []
    g = []
    c = []

    for cluid in backbone_node[1:]:
        if cluid == ".":
            size = 0
            gene = "."
            chromosome = "."
        else:
            query_clu = cluid.split(",")
            size_p = []
            gene_p = []
            chromosome_p = []

            for subcluid in query_clu:
                try:
                    query = clu_data[ind][int(subcluid)-1].split()
                    size_p.append(int(query[3]))
                    gene_p.append(query[2])
                    chromosome_p.append(query[4])
                except:
                    size_p.append(-9999)
                    gene_p.append(str(subcluid)+":?")
                    chromosome_p.append(str(subcluid)+":?")

            split_str = ','
            size = sum(size_p)
            gene = split_str.join(list(set(gene_p)))
            chromosome = split_str.join(list(set(chromosome_p)))

            size_csv_file.write(files[ind+1] + "\t" + backbone_node[0] + "\t" + str(size) + "\n")

        ind += 1
        s.append(str(size))
        g.append(str(gene))
        c.append(str(chromosome))

    Chr = '|'.join(list(set(c)))
    try:
        Clu = backbone_node[0]
    except IndexError:
        size_csv_file.close()
        size_file.close()
        clu_file.close()
        chr_file.close()
        chr_sum.close()
        sys.exit()

    chr_sum.write(Clu + '\t' + Chr + '\n')
    subsize = Clu + '\t' + "\t".join(s) + "\n"
    size_file.write(subsize)
    subclu = Clu + '\t' + "\t".join(g) + "\n"
    clu_file.write(subclu)
    subchr = Clu + '\t' + "\t".join(c) + "\n"
    chr_file.write(subchr)

size_csv_file.close()
size_file.close()
clu_file.close()
chr_file.close()
chr_sum.close()
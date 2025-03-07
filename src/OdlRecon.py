#!/usr/bin/env python3
#coding=utf-8
import argparse
import numpy as np
import pandas as pd
import re

# Add command line arguments
parser = argparse.ArgumentParser(description='Infer ODL events based on dynamic programming', add_help=False)
parser.add_argument('-h', '--help', action='store_true', help='Show this help message and exit')
parser.add_argument('-p', '--path', required=True, help='Path to input/output files')
parser.add_argument('-g', '--genome_info', required=True, help='Path to genome info file')
parser.add_argument('-s', '--hoc_sizes', required=True, help='Path to hoc sizes file')
parser.add_argument('-t', '--speices_tree', required=True, help='Path to speices tree file')
parser.add_argument('-o', '--output_prefix', required=True, help='Prefix for output files')

args = parser.parse_args()

# Read leaf node file (Newick format)
with open(args.speices_tree, 'r') as f:
    tree_str = f.read().strip()
    # Extract leaf node names from Newick format tree
    selected_groups = re.findall(r'([A-Za-z0-9_]+):', tree_str)[::-1]

## Read genome information from file
genome_info = pd.read_csv(args.genome_info, sep='\t', header=None, names=['spec', 'group', 'genome'])
group_dict = genome_info[genome_info['group'].isin(selected_groups)].set_index('genome')['group'].to_dict()

# Read genome size information
genome_sizes = pd.read_csv(args.hoc_sizes, sep='\t', header=0)
hocs = genome_sizes['hoc'].unique()

# Create output files
median_file = open(f"{args.path}/{args.output_prefix}.median.tsv", 'w')
node_odl_file = open(f"{args.path}/{args.output_prefix}.node.odl.tsv", 'w')
leaf_odl_file = open(f"{args.path}/{args.output_prefix}.leaf.odl.tsv", 'w')

for hoc, lines in genome_sizes.groupby('hoc'):
    not_zero_median_dict = {}
    for _, row in lines.iterrows():
        genome = row['genome']

        if genome in group_dict:
            hoc_size = row['gene_count']
            if int(hoc_size)!= 0:
                if group_dict[genome] not in not_zero_median_dict:
                    not_zero_median_dict[group_dict[genome]] = {genome: int(hoc_size)}
                else:
                    not_zero_median_dict[group_dict[genome]][genome] = int(hoc_size)

    for group in not_zero_median_dict:
        not_zero_median_dict[group] = sorted(not_zero_median_dict[group].items(), key=lambda x: x[1])[len(not_zero_median_dict[group])//2]

    median_file.write(hoc + '\t' + str(not_zero_median_dict) + '\n')

    leaf_hoc_size = []
    for group in selected_groups:
        if group in not_zero_median_dict:
            leaf_hoc_size.append(not_zero_median_dict[group][1])
        else:
            leaf_hoc_size.append(0)

    # Solve using dynamic programming
    def dynamic_programming_solve(L):
        n = len(L) - 1  # 10
        dp = np.full((n+1, n+1), float('inf'))
        path = np.zeros((n+1, n+1), dtype=int)

        # Initialization
        for i in range(n+1):
            dp[0, i] = abs(L[0] - L[i])

        # Fill the dynamic programming table
        for i in range(1, n):
            for j in range(i, n+1):
                for k in range(i-1, j):
                    cost = dp[i-1, k] + abs(L[i] - L[j]) + abs(L[k] - L[j])
                    if cost < dp[i, j]:
                        dp[i, j] = cost
                        path[i, j] = k

        # Backtrack to find the optimal solution
        X = [0] * n
        X[n-1] = n
        for i in range(n-1, 0, -1):
            X[i-1] = path[i, X[i]]

        # Calculate the final X values
        for i in range(n):
            X[i] = L[X[i]]

        return X

    L = leaf_hoc_size[::-1]  # Reverse order
    X = dynamic_programming_solve(L)

    tag_dict = {}
    if X[0] > 1:
        tag_dict['sp0'] = 'origin|duplicate'
    elif X[0] > 0:
        tag_dict['sp0'] = 'origin'
    else:
        tag_dict['sp0'] = '.'

    for i in range(len(X) - 1):
        diff = X[i + 1] - X[i]
        if diff > 1 and X[i] == 0:
            tag_dict['sp' + str(i + 1)] = 'origin|duplicate'
        elif diff > 0 and X[i] == 0:
            tag_dict['sp' + str(i + 1)] = 'origin'
        elif diff > 0:
            tag_dict['sp' + str(i + 1)] = 'duplicate'
        elif diff < 0 and X[i + 1] == 0:
            tag_dict['sp' + str(i + 1)] = 'loss'
        elif diff < 0:
            tag_dict['sp' + str(i + 1)] = 'shrink'
        else:
            tag_dict['sp' + str(i + 1)] = '.'

    tag_dict2 = {}
    for i in range(len(X)):
        diff = L[i] - X[i]
        if diff > 0:
            tag_dict2['L' + str(i + 1)] = 'duplicate'
        elif diff < 0 and L[i] == 0:
            tag_dict2['L' + str(i + 1)] = 'loss'
        elif diff < 0:
            tag_dict2['L' + str(i + 1)] = 'shrink'
        else:
            tag_dict2['L' + str(i + 1)] = '.'

    # Handle special case for the last element
    diff = L[9] - X[8]
    if diff > 0:
        tag_dict2['L10'] = 'duplicate'
    elif diff < 0 and L[9] == 0:
        tag_dict2['L10'] = 'loss'
    elif diff < 0:
        tag_dict2['L10'] = 'shrink'
    else:
        tag_dict2['L10'] = '.'
    
    node_odl_file.write(hoc + '\t' + '\t'.join(tag_dict.values()) + '\n')
    leaf_odl_file.write(hoc + '\t' + '\t'.join(tag_dict2.values()) + '\n')

median_file.close()
node_odl_file.close()
leaf_odl_file.close()
#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-03-08, 15:07:33
# modified date: 2024-03-08, 15:07:33

'''
Jinxin Meng, 2024
filter profile by coverage in another file.
Required files:
1. reads count profile
2. breadth of coverage profile
'''

import sys

if len(sys.argv) != 5:
    print("Usage: profile_filter_by_cvg.py [*.rc] [*.cvg] [breath_cutoff] [out_file]")
    sys.exit()

out_f = open(str(sys.argv[4]), "w")
with open(str(sys.argv[1]), "r") as f1, open(str(sys.argv[2]), "r") as f2:
    header = f1.readline()
    out_f.write(header)
    next(f2)
    for l1, l2 in zip(f1, f2):
        l1 = l1.strip().split()
        l2 = l2.strip().split()
        name = l1[0]
        x = [int(x) for x in l1[1:]]
        y = [int(float(x) > int(sys.argv[3])) for x in l2[1:]]
        res = [i * j for i, j in zip(x, y)]
        if sum(res) != 0:
            res = [name] + [str(x) for x in res]
            out_f.write('\t'.join(res) + '\n')
out_f.close()

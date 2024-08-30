#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-03-08, 15:35:09
# modified date: 2024-03-08, 17:22:26

'''
Jinxin Meng, 2024
summary profile by the specified group.
Required files: rc|tpm|cvg|.. profile
default group on first column, and method is sum.
'''

import sys, re

if len(sys.argv) != 3:
    print("Usage: profile_summary_by_group.py [*.rc|tpm|cvg ..] [out_file]")
    sys.exit()

out_f = open(str(sys.argv[2]), "w")
with open(str(sys.argv[1]), "r") as f1:
    header = f1.readline()
    out_f.write(header)
    res_dict = {}
    for l in f1:
        l = l.strip().split()
        name = re.findall("(PUL\d+)", l[0])[0]
        if name not in res_dict:
            res_dict[name] = [int(x) for x in l[1:]]
        else:
            res_dict[name] = [x + int(y) for x, y in zip(res_dict[name], l[1:])]
for k, v in res_dict.items():
    out_f.write("\t".join([k] + [str(x) for x in v]) + "\n")
out_f.close()

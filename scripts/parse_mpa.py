#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-05-17, 22:16:50
# modified date: 2024-05-17, 22:16:50

import sys, argparse, re

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", metavar="in_f", type=str, help="metaphlan profile")
    parser.add_argument("-o", metavar="out_f", type=str, help="out prefix")
    parser.add_argument("-l", metavar="level", type=str, default="p", help="specify one or multiple (comma-delimited characters) taxonomic level to output [p(default)|c|o|f|g|s|t]")
    parser.add_argument("--short", action="store_true", help="profile name only retain the current level")
    parser.add_argument("--out_all", action="store_true", help="output all levels profile")
    parser.add_argument("--taxa", action="store_true", help="force to output taxonomy file")
    parser.add_argument("--split", action="store_true", help="split taxonomy infomation")
    parser.add_argument("--add_header", action="store_true", help="add header to taxonomy file")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args

def store_query():
    query = {"k": "kingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species", "t": "strain"}
    return query

def determine_taxa(taxa_list):
    query = store_query()
    subject = ", ".join(taxa_list)
    retain = []
    for i in query.keys():
        x = i + "__"
        if re.search(x, subject):
            retain.append(i)
    retain = [ query[x] for x in query.keys()]
    return retain

def parse_taxa(profile, out_f, split, add_header):
    init_tax = []
    with open(profile, "r") as f:
        for i in f:
            x = i.strip().split("\t")[0]
            if re.search("name", x):
                continue
            init_tax.append(x)
    init_subject = ",".join(init_tax)
    final_tax = []
    for i in init_tax:
        if len(re.findall(re.escape(i), init_subject)) == 1:
            final_tax.append(i)
    out_f = open(out_f + ".taxa", "w")
    if split is True:
        if add_header is True:
            out_f.write("\t".join(determine_taxa(final_tax)) + "\n")
        for i in final_tax:
            out_f.write("\t".join(i.split("|")) + "\n")
    else:
        for i in final_tax:
            out_f.write(i + "\n")
    out_f.close()

def parse_profile(profile, level, out_f, out_short):
    out_f = open(out_f, "w")
    pattern = re.compile("(" + level + "__\w+(?!\|))$")
    with open(profile, "r") as f:
        out_f.write(f.readline())
        for i in f:
            l = i.strip().split("\t", 1)
            if re.search(pattern, l[0]):
                if out_short is True:
                    name = re.findall(pattern, l[0])[0]
                    out_f.write(name + "\t" + l[1] + "\n")
                else:
                    out_f.write(i)
    out_f.close()

def main(in_f, out_f, level, out_short, out_all, out_taxa, split, add_header):
    if out_taxa is True:
        parse_taxa(in_f, out_f, split, add_header)
        sys.exit(0)
        
    query = store_query()
    if out_all is True:
        level = "p,c,o,f,g,s"
    for i in level.split(","):
        if i not in query.keys():
            sys.exit("The parameter error in -l.")
        parse_profile(in_f, i, out_f + "." + query[i] + ".tsv", out_short)

if __name__ == "__main__":
    args = get_args()
    main(args.i, args.o, args.l, args.short, args.out_all, args.taxa, args.split, args.add_header)


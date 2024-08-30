#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-04-02, 22:54:53
# modified date: 2024-04-03, 00:13:01

import argparse, sys
import pandas as pd
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu # Mann-Whitney U test (sometimes called the Wilcoxon rank-sum test)

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-profile", metavar="profile", type=str, default="", help="in_f profile-like data.")
    parser.add_argument("-metadata", metavar="metadata", type=str, default="", help="in_f sample-group data.")
    parser.add_argument("-group_pair", metavar="group_pair", type=str, default="", help="a group pair file contain two field with tab-delimited (x|y)")
    parser.add_argument("-center_group", metavar="center_group", type=str, default="", help="a group name used to comared with all another group")
    parser.add_argument("-method", metavar="method", type=str, default="wilcox", help="method in [t, wilcox]. default: wilcox")
    parser.add_argument('-out', metavar="out_f", type=str, default="", help='output file')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args

def print_process_bar(iteration, total, prefix="", suffix="", decimals=1, length=50, fill="█"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    fill_length = int(length * iteration // total)
    bar = fill * fill_length + "-" * (length - fill_length)
    sys.stdout.write('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
    if iteration == total:
        sys.stdout.write("\n")
        sys.stdout.flush()

def get_gp(group, center_group):
    group_list = list(set(group.loc[:, "group"]))
    out_list = []
    if center_group == "":
        for i in range(len(group_list[:-1])):
            for j in range(i + 1, len(group_list)):
                out_list.append([group_list[i], group_list[j]])
    elif center_group != "" and center_group in group_list:
        group_list.remove(center_group)
        for i in group_list:
            out_list.append([i, center_group])
    else:
        sys.exit("center_group not in group_list ..")
    return(out_list)

def parse_gp(group_pair):
    gp = []
    with open(group_pair, "r") as f:
        for i in f:
            l = i.strip().split("\t")
            gp.append([l[0], l[1]])
    return gp

def convert_plab(pval):
    if pval < 0.001:
        plab = "***"
    elif pval < 0.01:
        plab = "**"
    elif pval < 0.05:
        plab = "*"
    else:
        plab = ""
    return plab

def main(profile, group, group_pair, center_group, method, out):
    profile = pd.read_csv(profile, sep="\t", header=0, index_col=0)
    group = pd.read_csv(group, sep="\t", header=0, index_col=None)

    if group_pair == "":
        gp = get_gp(group, center_group)
    else:
        gp = parse_gp(group_pair)

    name = profile.index
    total = len(name) * len(gp)
    iteration = 1
    if method == "t":
        print_process_bar(0, total, prefix="Processing:", suffix="Complete", length=50)
        out_f = open(out, "w")
        out_f.write("name\tgroup_pair\tpval\tplab\tmethod\n")
        for i in name:
            for j in gp:
                sample_x = group.loc[group["group"]==j[0], "sample"]
                sample_y = group.loc[group["group"]==j[1], "sample"]
                vec_x = profile.loc[i, sample_x]
                vec_y = profile.loc[i, sample_y]
                t_statistic, p_value = ttest_ind(vec_x, vec_y)
                plab = convert_plab(p_value)
                out_f.write("\t".join([i, j[0] + "_vs_" + j[1], str(round(p_value, 4)), plab, method]) + "\n")
                print_process_bar(iteration, total, prefix="Processing:", suffix="Complete", length=50)
                iteration += 1
        out_f.close()
    elif method == "wilcox":
        print_process_bar(0, total, prefix="Processing:", suffix="Complete", length=50)
        out_f = open(out, "w")
        out_f.write("name\tgroup_pair\tpval\tplab\tmethod\n")
        for i in name:
            for j in gp:
                sample_x = group.loc[group["group"]==j[0], "sample"]
                sample_y = group.loc[group["group"]==j[1], "sample"]
                vec_x = profile.loc[i, sample_x]
                vec_y = profile.loc[i, sample_y]
                w_statistic, p_value = mannwhitneyu(vec_x, vec_y)
                plab = convert_plab(p_value)
                out_f.write("\t".join([i, j[0] + "_vs_" + j[1], str(round(p_value, 4)), plab, method]) + "\n")
                print_process_bar(iteration, total, prefix="Processing:", suffix="Complete", length=50)
                iteration += 1
        out_f.close()
    else:
        sys.exit("method should specify t or wilcox ..")

if __name__ == "__main__":
    args = get_args()
    main(args.profile, args.metadata, args.group_pair, args.center_group, args.method, args.out)


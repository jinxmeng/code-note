#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-03-01, 18:15:11
# modified date: 2024-03-13, 19:37:35

# 2024-03-13: add --no_header parameter.

import argparse, sys, signal 
from io import StringIO
import pandas as pd

signal.signal(signal.SIGPIPE, signal.SIG_DFL) # 忽略 SIGPIPE 信号

def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", metavar="in_f", type=str, default="", help="input table or stdin")
    parser.add_argument("-s", metavar="sep", type=str, default="\t", help="what separate for table")
    parser.add_argument("-n", metavar="name", type=str, default="", help="only column or row of this name to stat")
    parser.add_argument("-f", metavar="field", type=int, default=-1, help="only column or row of this index to stat")
    parser.add_argument("-m", metavar="method", type=str, default="all", help="method in [n, sum, mean, std, min, q_25, q_50, q_75, max, all]. default: all")
    parser.add_argument("-r", action="store_true", help="row for profile stat. default: column")
    parser.add_argument('-o', metavar="out_f", type=str, default="", help='output file or stdout')
    parser.add_argument('--no_header', action="store_true", help='in_f without header')
    args = parser.parse_args()
    return args

def f_stat(df, method):
    stat = {"name": df.columns}
    if method == "n":
        stat["n"] = df.count()
    elif method == "sum":
        stat["sum"] = df.sum()
    elif method == "mean":
        stat["mean"] = round(df.mean(), 4)
    elif method == "std":
        stat["std"] = round(df.std(), 4)
    elif method == "min":
        stat["min"] = df.min()
    elif method == "q_25":
        stat["q_25"] = df.quantile(0.25)
    elif method == "q_50":
        stat["q_50"] = df.quantile(0.5)
    elif method == "q_75":
        stat["q_75"] = df.quantile(0.75)
    elif method == "max":
        stat["max"] = df.max()
    else:
        stat["n"] = df.count()
        stat["sum"] = df.sum()
        stat["mean"] = round(df.mean(), 4)
        stat["std"] = round(df.std(), 4)
        stat["min"] = df.min()
        stat["q_25"] = df.quantile(0.25)
        stat["q_50"] = df.quantile(0.5)
        stat["q_75"] = df.quantile(0.75)
        stat["max"] = df.max()
    stat = pd.DataFrame(stat)
    return stat

def s_stat(series, method):
    stat = {"name": series.name}
    if method == "n":
        stat["n"] = series.count()
    elif method == "sum":
        stat["sum"] = series.sum()
    elif method == "mean":
        stat["mean"] = round(series.mean(), 4)
    elif method == "std":
        stat["std"] = round(series.std(), 4)
    elif method == "min":
        stat["min"] = series.min()
    elif method == "q_25":
        stat["q_25"] = series.quantile(0.25)
    elif method == "q_50":
        stat["q_50"] = series.quantile(0.5)
    elif method == "q_75":
        stat["q_75"] = series.quantile(0.75)
    elif method == "max":
        stat["max"] = series.max()
    else:
        stat["n"] = series.count()
        stat["sum"] = series.sum()
        stat["mean"] = round(series.mean(), 4)
        stat["std"] = round(series.std(), 4)
        stat["min"] = series.min()
        stat["q_25"] = series.quantile(0.25)
        stat["q_50"] = series.quantile(0.5)
        stat["q_75"] = series.quantile(0.75)
        stat["max"] = series.max()
    stat = pd.DataFrame(stat, index=[0])
    return stat

def main(in_f, sep, name, field, method, row, out_f, no_header):
    if no_header:
        header = None
    else:
        header = 0
    if in_f != "":
            df = pd.read_csv(in_f, sep=sep, header=header, index_col=0, skiprows=0)
    else:
        print("Reading from stdin ..", file=sys.stderr)
        df = sys.stdin.read()
        if df:
            df = StringIO(df)
            df = pd.read_csv(df, sep=sep, header=header, index_col=0, skiprows=0)
        else:
            print("can not receive data ..", file=sys.stderr)
    if row:
        df = df.T
    else:
        df = df
    if name != "":
        series = df[name]
        res = s_stat(series, method=method)
    if field != -1:
        series = df.iloc[:, field - 1]
        res = s_stat(series, method=method)
    if name == "" and field == -1:
        res = f_stat(df, method=method)
    if out_f != "":
        res.to_csv(out_f, sep="\t", index=False, header=True)
    else:
        print(res.to_string(index=False), file=sys.stdout)

if __name__ == '__main__':
    args = get_args()
    main(args.i, args.s, args.n, args.f, args.m, args.r, args.o, args.no_header)

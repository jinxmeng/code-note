#!/share/software/anaconda3/bin/python
# -*- encoding: utf-8 -*-
# Created  date :  2023-03-01
# Modiffed date :  2021-03-01

import gc
import json
import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', metavar='in', type=str, help='Input file')
    parser.add_argument('-o', metavar='out', type=str, help='Output file')
    parser.add_argument('-w', metavar='key', type=str, help='Specify key words. Example: mRNA; mRNA,gene ..')
    args = parser.parse_args()
    return args

def getdictvalue(d, co_de):
    result = []
    if isinstance(d, dict):
        try:
            value = d[co_de]
            result.append(value)
        except Exception as e:
            pass
        for valuedd in d.values():
            if isinstance(valuedd, dict):
                yied_result = getdictvalue(valuedd, co_de)
                if len(yied_result) != 0:
                    result.append(getdictvalue(valuedd, co_de))
            elif isinstance(valuedd, (list, tuple)):
                for item in d:
                    valueitem = getdictvalue(item, co_de)
                    if valueitem != "None" and valueitem is not None and len(valueitem) != 0:
                        if valueitem not in result:
                            result.append(valueitem)
            elif isinstance(valuedd,str):
                try:
                    resultd = eval(valuedd)
                    for item in resultd:
                        value = getdictvalue(item, co_de)
                        if value != "None" and value is not None and len(value) != 0:
                            if value not in result:
                                result.append(value)
                except:
                    pass
        gc.collect()

    elif isinstance(d, (list, tuple)):
        for item in d:
            value = getdictvalue(item, co_de)
            if value != "None" and value is not None and len(value) != 0:
                if value not in result:
                    result.append(value)
        gc.collect()
    else:
        try:
            d=eval(d)
            for item in d:
                value = getdictvalue(item, co_de)
                if value != "None" and value is not None and len(value) != 0:
                    if value not in result:
                        result.append(value)
            gc.collect()
        except:
            pass
    return result

def unlist(complist):
    while True:
        if isinstance(complist, list):
            complist = complist[0]
        else:
            break
    return complist

if __name__ == "__main__":
        args = get_args()
        words = args.w.split(sep = ',')
        f = open(args.i,  "r", encoding='utf-8')
        f2 = open(args.o, "a", encoding='utf-8') 
        f2.write('\t'.join(words))
        f2.write('\n')
                    
        lines = f.readlines()
        for i in lines:
            dat = json.loads(i) # json.load执行一个文件， json.loads执行文件的每一行
        
            # 新建一个list存储数据                          
            lst = []
            for j in words:
                lst.append(str(unlist(getdictvalue(dat, j))))
            f2.write('\t'.join(lst))
            f2.write('\n')
        f.close()

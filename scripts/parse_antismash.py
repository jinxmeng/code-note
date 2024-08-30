#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-01-23, 23:20:13
# modified date: 2024-05-17, 10:49:42

import re, sys, argparse
from bs4 import BeautifulSoup

if len(sys.argv) != 3:
    sys.exit("Usage: python3 " + sys.argv[0] + " [in_f index.html] [out_prefix *tsv]")

def get_contigs(soup):
    info = soup.find_all('div', attrs={'class':'record-overview-header'})
    contigs = []
    for i in info:
        for j in i.find_all('strong'):
            x = j.text.strip()
            contigs.append(x)
    return contigs

def get_res(soup, in_f):
    info = soup.find_all('table', attrs={'class':'region-table'})[0:-1]  # 根据属性来查找元素使用attr参数
    if len(info) == 0:
        sys.exit('Not find BGCs in ' + in_f + '!')
    header = []
    for i in info[0].find_all('th'):
        header.append(i.text)
    res = []
    count = []
    for i in info:
        x = i.find_all('td')
        row = []
        n = 1
        for j in x:
            y = j.text.strip()
            if y == "":
                continue
            if re.search(r'&nbsp', y):# re.search 判断是否字符串中有某个或者某些字符
                y = y.replace('&nbsp', ' ')
            if re.search('Region', y) and len(row) != 0:
                res.append(row)
                row = []
                n += 1
            row.append(y)
        res.append(row)
        count.append(n)
    res2 = []
    for i in res:
        if len(i) > 4:
            x = i[0:4]
            x.append(' - '.join(i[4:-1]))
            x.append(i[-1])
            res2.append(x)
        else:
            res2.append(i)
    return header, res2, count

def main(in_f, out_f):
    html = open(in_f, 'r')
    soup = BeautifulSoup(html, 'html.parser')
    header, res, count = get_res(soup, in_f)
    contigs = get_contigs(soup)
    contigs2 = []
    for i in range(len(contigs)):
        contigs2 += [contigs[i]] * count[i]
    out_f = open(out_f, 'w')
    out_f.write('contigs' + '\t' + '\t'.join(header) + '\n')
    for i in range(len(res)):
        out_f.write(contigs2[i] + '\t' + '\t'.join(res[i]) + '\n')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

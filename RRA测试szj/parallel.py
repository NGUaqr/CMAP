import subprocess
import sys
import os
import PIL

from multiprocessing import Pool

# 调用cmd命令的例子
Lfile=sys.argv[1]#lincs drug profile
Sfile=sys.argv[2]#disease signature

def f(name):
    if name==1:
        os.system('python SS_fisher.py '+Lfile+' '+Sfile)
    if name==2:
        os.system('python SS_corr.py '+Lfile+' '+Sfile)
    if name==3:
        os.system('python SS_LINCS.py '+Lfile+' '+Sfile)
    if name==4:
        os.system('python SS_ZhangScore.py '+Lfile+' '+Sfile)
    if name==5:
        os.system('perl SS_CMAP.pl '+Lfile+' '+Sfile+' 0.30')
if __name__ == '__main__':
    with Pool(5) as p:
        p.map(f, [1,2,3,4,5])


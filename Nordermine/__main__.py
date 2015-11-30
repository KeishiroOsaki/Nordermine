'''
Created on 2015/11/30

@author: osakikeishiro
'''

import csv
import sys

if __name__ == '__main__':
    print("おはようございます おはよう世界")
    args = sys.argv
    k = int(args[1])
    inputCSVPath = args[2]
    outputCSVPath = args[3]
    useColumns = []
    for i in range(4,len(args)):
        useColumns.append(int(args[i]))
    print(useColumns)
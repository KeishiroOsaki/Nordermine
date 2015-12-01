# -*- coding: utf-8 -*-
'''
Created on 2015/11/30

@author: osakikeishiro

required Python3.5
'''

import csv
import sys
from collections import Counter
from scipy import special
import math
import random
import time

def logsumexp(x, y, flg):
    if flg != 0:
        return y
    if x == y:
        return x + 0.69314718055
    vmin = 0
    vmax = 0
    if x > y:
        vmax = x
        vmin = y
    else:
        vmax = y
        vmin = x
    if vmax > vmin + 50:
        return vmax
    else:
        return vmax + math.log(math.exp(vmin - vmax) + 1.0)

class Nordermine:
    def __init__(self, records, k, useColumns,header):
        self.records = records
        self.M = len(useColumns)
        self.z = [-1 for i in range(0, len(records))]
        self.k = k
        self.useColumns = useColumns
        self.header = header
        
        # 重複を除いた値のリストを作る
        self.dValues = []
        self.numofDvalues = []
        for i in range(0, self.M):
            tmplist = []
            for j in range(0, len(rows)):
                tmplist.append(rows[j][i])
            self.dValues.append(sorted(set(tmplist)))
            self.numofDvalues.append(len(self.dValues[i]))
        
        # 値からインデックスを求めるための辞書を生成
        self.dValuesDict = []
        for i in range(self.M):
            tmpDict = {}
            for j in range(len(self.dValues[i])):
                tmpDict[self.dValues[i][j]] = j
            self.dValuesDict.append(tmpDict)
        
        #レコードをインデックス化
        self.indexRecords = []
        for i in range(len(self.records)):
            tmp = []
            for j in range(self.M):
                tmp.append(self.dValuesDict[j][self.records[i][j]])
            self.indexRecords.append(tmp)
        
        # ハイパーパラメータの初期化
        self.hyperParameters = [0.1 for j in range(0, self.M)]
        self.hyperParameters[0] = 50.0 / self.k
        
        self.Nk = [0 for j in range(0, self.k)]
        self.Nu = [0 for j in range(0, len(self.dValues[0]))]
        
        tmp = [self.records[i][0] for i in range(len(self.records))]
        counter = Counter(tmp)
        for item, cnt in counter.most_common():
            self.Nu[self.dValuesDict[0][item]] = cnt
        
        # エントリの割当回数を記憶するための行列
        self.Nxk = []
        for i in range(0, self.M):
            self.Nxk.append([[0 for i in range(self.k)] for j in range(len(self.dValues[i]))])
        
        self.ansModel = []
        for i in range(0, self.M):
            self.ansModel.append([[0.0 for i in range(self.k)] for j in range(len(self.dValues[i]))])
   
    def updateParameter(self):
        eps = 1e-8
        nume = -self.numofDvalues[0] * self.k * special.digamma(self.hyperParameters[0])
        deno = -self.numofDvalues[0] * self.k * special.digamma(self.hyperParameters[0] * self.k)
        for i in range(self.numofDvalues[0]):
            deno += self.k * special.digamma(self.Nu[i] + self.hyperParameters[0] * self.k)
            for j in range(self.k):
                nume += special.digamma(self.Nxk[0][i][j] + self.hyperParameters[0])
        self.hyperParameters[0] *= nume / deno
        if self.hyperParameters[0] < eps:
            self.hyperParameters[0] = eps
            
        for g in range(1, len(self.hyperParameters)):
            nume = -self.numofDvalues[g] * self.k * special.digamma(self.hyperParameters[g])
            deno = -self.numofDvalues[g] * self.k * special.digamma(self.hyperParameters[g] * self.numofDvalues[g])
            for i in range(self.k):
                deno += self.numofDvalues[g] * special.digamma(self.Nk[i] + self.hyperParameters[g] * self.numofDvalues[g])
                for j in range(self.numofDvalues[g]):
                    nume += special.digamma(self.Nxk[g][j][i] + self.hyperParameters[g])
            self.hyperParameters[g] *= nume / deno
            if self.hyperParameters[g] < eps:
                self.hyperParameters[g] = eps

    def sampling(self,index):
        posts = []
        
        rowIdx = self.indexRecords[index]
        
        topic = self.z[index]
        if topic != -1: #2回目以降のiterの場合に実行
            self.Nk[topic] -= 1
            for i in range(self.M):
                self.Nxk[i][rowIdx[i]][topic] -= 1
            
        for l in range(self.k):
            tmp = 1.0
            
            nume=self.Nxk[0][rowIdx[0]][l] + self.hyperParameters[0]
            deno = self.Nu[rowIdx[0]] + self.hyperParameters[0] * self.k
            tmp *= nume/deno
            
            for m in range(1,self.M):
                nume = self.Nxk[m][rowIdx[m]][l] + self.hyperParameters[m]
                deno = self.Nk[l] + self.hyperParameters[m] * len(self.dValues[m])
                tmp *= nume/deno
            posts.append(tmp)
            
        Z = sum(posts)
        u = Z * random.random()
        nextz = 0
        currprob = posts[0]
        
        while u > currprob:
            nextz += 1
            currprob += posts[nextz]
        
        self.z[index] = nextz
        self.Nk[nextz] += 1
        for m in range(self.M):
            self.Nxk[m][rowIdx[m]][nextz] += 1
        

    def inference(self,flg):
        maxStep = 5
        
        for i in range(maxStep):
            print("iter: " + str(i+1))
            for j in range(len(self.records)):
                self.sampling(j)
            print(self.Nk)
            if flg == True:
                self.updateParameter()
    def calcParams(self):
        for m in range(len(self.ansModel)):
            for k in range(self.k):
                for i in range(len(self.dValues[m])):
                    if m == 0:
                        self.ansModel[m][i][k] = (self.Nxk[m][i][k] + self.hyperParameters[m])/(self.Nu[i]+self.hyperParameters[m]*self.k)
                    else:
                        self.ansModel[m][i][k] = (self.Nxk[m][i][k] + self.hyperParameters[m])/(self.Nk[k]+self.hyperParameters[m]*len(self.dValues[m]))
    
    def writeParams(self,directory):
        for m in range(self.M):
            fileName = ""
            if m == 0:
                fileName = directory + "O"
            elif m == self.M-1:
                fileName = directory + "C"
            else:
                fileName = directory + "A" + str(m)
            f = open(fileName, 'w',encoding = 'UTF-8')
            writer = csv.writer(f, lineterminator='\n',delimiter=' ')
            writer.writerows(self.ansModel[m])
            f.close()
            
            f = open(fileName + ".tsv", 'w',encoding = 'UTF-8')
            writer = csv.writer(f, lineterminator='\n',delimiter='\t', quoting=csv.QUOTE_NONNUMERIC)
            hd = []
            hd.append(self.header[m])
            for i in range(self.k):
                hd.append("topic" + str(i +1))
            writer.writerow(hd)
            for i in range(len(self.ansModel[m])):
                r = []
                r.append(self.dValues[m][i])
                for j in range(self.k):
                    r.append(self.ansModel[m][i][j])
                writer.writerow(r)
            f.close()

if __name__ == '__main__':
    args = sys.argv
    k = int(args[1])
    inputCSVPath = args[2]
    outputPath = args[3]
    useColumns = []
    for i in range(4, len(args)):
        useColumns.append(int(args[i]) - 1)
    print(useColumns)
    
    print("データセット読み込み開始")
    rows = []
    f = open(inputCSVPath, 'r', encoding='utf-8')
    reader = csv.reader(f)
    header = next(reader)  # ヘッダーの読み飛ばし
    for row in reader:
        line = []
        for i in range(0, len(useColumns)):
            line.append(row[useColumns[i]])
        rows.append(line)
    print(header)
    print(rows[0])
    useHeader = []
    for i in range(len(useColumns)):
        print(header[useColumns[i]] , end = " ")
        useHeader.append(header[useColumns[i]])
    print("")
    print("読み込み完了 : " + str(len(rows)) + "件")
    
    model = Nordermine(rows, k, useColumns,useHeader)
    
    print("サンプリング開始")
    start = time.time()
    model.inference(True)
    elapsed_time = time.time() - start
    print (("elapsed_time:{0}".format(elapsed_time)) + "[sec]")
    model.calcParams()
    model.writeParams(outputPath)

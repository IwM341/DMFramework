import matplotlib.pyplot as plt
import numpy as np
import os
import sys

xlabel = ""
ylabel = ""
title = ""
res = 100
isLegend = False
yscale = "linear"
Arrays = []
fnames = []

i = 1
N = len(sys.argv)
while(i < N ):
    if(sys.argv[i] == "-xl"):
        if(i+1 < N):
            i+=1
            xlabel = sys.argv[i]
        else:
            print("expect xlabel after -xl")
        i+=1
    elif(sys.argv[i] == "-yl"):
        if(i+1 < N):
            i+=1
            ylabel = sys.argv[i]
        else:
            print("expect ylabel after -yl")
        i+=1
    elif(sys.argv[i] == "-tl"):
        if(i+1 < N):
            i+=1
            title = sys.argv[i]
        else:
            print("expect title after -tl")
        i+=1
    elif(sys.argv[i] == "-res"):
        if(i+1 < N):
            i+=1
            try:
                res = int(sys.argv[i])
            except Exception:
                print("expect resolution number instead of " + sys.argv[i])
                print(Exception)
        else:
            print("expect ylabel after -yl")
        i+=1
    elif(sys.argv[i] == "-sl"):
        isLegend = True
        i+=1
    elif(sys.argv[i] == "log"):
        yscale = "log"
        i+=1
    
    else:
        filename = sys.argv[i]
        try:
            arr = np.loadtxt(filename).transpose()
            Arrays.append(arr)
            fnames.append(filename)
        except:
            print("couldnt load filename: " + filename)
        i+=1


plt.rcParams['figure.dpi'] = res

for arr in Arrays:
    plt.plot(arr[0],arr[1])

if(isLegend):
    plt.legend(fnames)

plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.yscale(yscale)

plt.show()
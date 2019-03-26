import sys
import os.path
import numpy as np

param = sys.argv

if len(param) > 1 and param[1]== "-h":
    print('usage: argc[1] ')
    print('argc[1] = total number of samples')
    sys.exit()

if len(param) != 2:
    print('wrong usage.')
    print('usage: argc[1] ')
    print('argc[1] = total number of samples')
    sys.exit()
    
#Check File
Setnum=int(param[1])

np.zeros(Setnum)
DataTmp=[]
DataEne=[]
DataC=[]

for i in range(0, int(Setnum)):
    str_ = "./output/SS_rand"+str(i)+".dat"
    if os.path.isfile(str_) != True:
        print('The file ' + str_ + ' does not exist.')
        sys.exit()
        
    f =  open(str_, 'r')
    count=0
    for line in f:
        if count ==0 or count ==1:
            count +=1
            continue    
        data = line.split()
        
        if(i ==0):
            DataTmp.append(np.zeros(Setnum))
            DataEne.append(np.zeros(Setnum))
            DataC.append(np.zeros(Setnum))
            
        DataTmp[count-2][i]=1.0/float(data[0])
        DataEne[count-2][i]=float(data[1])
        DataC[count-2][i]=pow(float(data[0]),2)*( float(data[2])-pow(float(data[1]),2))
        count+=1
    f.close()

count=0
for datatmp in DataTmp:
    print datatmp.mean(), datatmp.std(), \
          DataEne[count].mean(), DataEne[count].std(),\
          DataC[count].mean(), DataC[count].std()\
    
    count+=1
    

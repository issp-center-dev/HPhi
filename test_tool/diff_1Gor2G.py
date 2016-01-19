import sys

#argv: 0=
input_list = sys.argv[1:]

try:
    strInputFolder = input_list[0]
    strTargetFolder = input_list[1]
    strFileType = input_list[2]
except:
    print('input error: inputFolder, outputFolder, strFileType.')

if strFileType=='2G':
    strFileName='/zvo_cisajscktalt.dat'
elif strFileType=='1G':
    strFileName='/zvo_cisajs.dat'
else:
    print ('strFileType must be 1G or 2G.')
    
strInputFilePath=strInputFolder + strFileName
strTargetFilePath=strTargetFolder + strFileName

print('Origin File is '+strInputFilePath)
print('Target File is '+strTargetFilePath)

try:
    readFileIn = open(strInputFilePath, 'r')
    readFileTarget = open(strTargetFilePath, 'r')
    idiffErr =0
    for lineIn in readFileIn:
        dataIn=lineIn.split()
        dataTarget=readFileTarget.readline()
        dataTarget=dataTarget.split()
        icnt =0
        ierr =0
        delta2 =0
        for data in dataTarget:
            delta2 =float(data)-float(dataIn[icnt])
            delta2 = delta2**2
            if delta2 > abs(float(data))**2*pow(10.0, -14) :
                print('Mismatch for the following values')
                print(dataIn)
                print(dataTarget)
                idiffErr +=1
            icnt +=1
    if idiffErr ==0:
        print('Files have same values!')
except:
    print('error for open files')

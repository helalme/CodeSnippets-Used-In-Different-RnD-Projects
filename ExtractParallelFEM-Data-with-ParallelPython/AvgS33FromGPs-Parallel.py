import sys
import os
import os.path
#import matplotlib.pyplot as plt
import pdb
from joblib import Parallel, delayed
import multiprocessing as mp
     

#num_cores = multiprocessing.cpu_count() 
def extractSig33(i,fpL):
    n=0                #Number of element
    #fileLinebyline=[]   #a text file will be saved as list of lines
    time=[]             #timesteps availale in the Ofiles
    gpStressAverage=[]  #averaged stress of a component based on all GPs from a single Ofile
    lineNumberTrack=0   #necessary to know expected data stays in which line, or after which line
    gpTracker=0         #Gauss point tracker, will be 0 at the beginning of each timestep
    gpStressSum=0       #Sum of a stress component from all GPs
    gpCount=0           #Counting gauss points
    #filePath=filePath[:-1]+ str(fp)
    with open(fpL)as f:
        #fileLineByLine=f.readlines()
        lineCount = -1
        for line in f:#ileLineByLine:
            lineCount += 1
            listFromLine=line.split()
            if len(listFromLine)>=4:
                #reading number of elements in an Ofile
                if listFromLine[0]=='Number' and listFromLine[1]=='of' and listFromLine[2]=='Elements':
                    n=int(listFromLine[12])
                if listFromLine[0]=='Computing' and listFromLine[1]=='solution':
                    strTime=listFromLine[4][:-1]
                    time.append(float(strTime))
                #calculating stress average from Gauss points
                if len(listFromLine)==5 and len(time)>=1 and listFromLine[0].isdigit()and listFromLine[1].isdigit():
                    #listFromLine[0] and listFromLine[1] are number of elements and materials number, respectively
                    if int(listFromLine[0])>=1 and int(listFromLine[0])<=n and int(listFromLine[1])>=1 and int(listFromLine[1])<=2:
                        gpTracker=lineCount
                if lineCount==gpTracker+1:
                    gpStressSum += float(listFromLine[2])
                    gpCount += 1
                if (gpCount==8*n) and (len(time)>=1) and (n>0):
                    #if (lineCount==gpTracker+4 or lineCount==gpTracker+20 or lineCount==gpTracker+26) and (listFromLine[0]=='*Command' or listFromLine[0]=='Rank['):
                    gpStressAverage.append(gpStressSum/(8*n))
                    gpStressSum=0
                    gpCount=0

    return [n,time,gpStressAverage]
                    
     

if __name__ == "__main__":
    cwd=os.getcwd()
    filePathList=[]
    filePath=cwd+'/O10grains_0001' #1st Ofile from Core-1

    #counting number of available Ofiles and allpaths saving
    OfileNos=0 #number of Ofiles
    while os.path.isfile(filePath):
        filePathList.append(filePath)
        OfileNos += 1
        if (OfileNos==9):
            filePath=filePath[:-1]
            filePath=filePath[:-1]+str(OfileNos+1)
        else:
            filePath=filePath[:-1]+str(OfileNos+1)

            
    pool=mp.Pool(4)
    '''
    embarrassingly parallel of data extracting from all Ofiles using Joblib, Parallel and delayed
    n_Timesteps_gpStressAverageFromAllOfiles = Parallel(n_jobs=4)(delayed(extractSig33)(i,filePathList[i-1] ) for i in range(1,OfileNos+1))
    '''
    #embarrassingly parallel of data extracting from all Ofiles using multiprocessing.Pool().apply() 
    n_Timesteps_gpStressAverageFromAllOfiles = [pool.apply(extractSig33, args=(i,filePathList[i-1])) for i in range(1,OfileNos+1)]
    pool.close()
            
    totalElements=0
    for i in range(OfileNos):
        totalElements += n_Timesteps_gpStressAverageFromAllOfiles[i][0]
        
    #pdb.set_trace()
    with open("StressAvgFromGPs-Parallel.txt","w") as of:
        sys.stdout=of
        for i in range(len(n_Timesteps_gpStressAverageFromAllOfiles[0][1])): 
            volStressSum=0
            for j in range(len(n_Timesteps_gpStressAverageFromAllOfiles)):
                volStressSum += n_Timesteps_gpStressAverageFromAllOfiles[j][2][i]*n_Timesteps_gpStressAverageFromAllOfiles[j][0]
            print(n_Timesteps_gpStressAverageFromAllOfiles[0][1][i], volStressSum/totalElements, sep='\t') 

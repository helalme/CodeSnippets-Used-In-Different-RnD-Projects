# At each timestep, extract desired data from different output files (Ofiles), e.g. stress and strain components along 33 directions
# Different Ofiles are coming from different cores from a parallel simulation
# Extracted data will be volume-averaged, and print into a file 

import sys
import os
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt



if __name__ == "__main__":
    cwd=os.getcwd()
    path=None
    Ofiles=None
    
    
    # User will run this script from command line, e.g. "pythons3 averageStressStrain.py O10grains_0001 8"
    if len(sys.argv)==1:
            print("InputError: Path for input data is missing.")
    if len(sys.argv)>1:
            path = cwd + "/" + sys.argv[1]
    if  len(sys.argv)>2:
            Ofiles = int(sys.argv[2])
    if len(sys.argv)>3:
            print("Argument should be 2 or 3")

    # Reading the complete file 
    stringList=[]
    if len(sys.argv)>=2:
        #path=cwd+'/O10grains'
        with open(path)as f:
            stringList=f.read().split()
        #print(1)

    '''
    # If the program is not run from command line
    path = cwd + '/O10grains_0001'
    Ofiles = 8
    with open(path)as f:
        stringList=f.read().split()
    ''' 	

    # Extracting Stress and strain components  at each time step     
    timeStrainStress=[]
    time=0
    for i in range(len(stringList)):
        if stringList[i]=='Computing' and stringList[i+1]=='solution':
            time=stringList[i+4][:-1] # [:-1] is for removing last char
        if stringList[i]=='Material' and stringList[i+1]=='All:':
            timeStrainStress.append([float(time),float(stringList[i+26]),float(stringList[i+14])])
    
    timestep=len(timeStrainStress)    
    if len(sys.argv)==3 or Ofiles is not None:
        for i in range(2,Ofiles+1):
            path=path[:-1]+ str(i)
            with open(path)as f:
                stringList=f.read().split()
            #print (i)
            for i in range(len(stringList)):
               if stringList[i]=='Computing' and stringList[i+1]=='solution':
                       time=stringList[i+4][:-1] # [:-1] is for removing last char
               if stringList[i]=='Material' and stringList[i+1]=='All:':
                       timeStrainStress.append([float(time),float(stringList[i+26]),float(stringList[i+14])])

    '''with open("strainStress.txt","w") as of:
        sys.stdout=of
        for i in timeStrainStress:
            print(*i, sep='\t')
	'''
    
    # Averaging stress and strain components and printing into file
    strainStressAvg=[]
    with open("avgStrainStress.txt","w") as of:
        sys.stdout=of
        for i in range(timestep):
            sumStrain=0
            sumStress=0
            for j in range(Ofiles):
                sumStrain=sumStrain+timeStrainStress[i+j*timestep][1]
                sumStress=sumStress+timeStrainStress[i+j*timestep][2]
            strainStressAvg.append([sumStrain/Ofiles,sumStress/Ofiles])
            print(timeStrainStress[i][0],sumStrain/Ofiles,sumStress/Ofiles,sep='\t') 
	
    # Plotting averaged stress and strain curve
    x,y=zip(*strainStressAvg)
    plt.plot(x,y)
    plt.xlabel('Average Strain (33-component)')
    plt.ylabel('Average Stress (33-component)')
    plt.show()       
   

    

    

   


  

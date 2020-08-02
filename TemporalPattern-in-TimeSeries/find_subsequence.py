import sys
import os
import pdb


#class containing two data members (attributes or properties) and and three function members (methods)  
class SubsequenceMaxSum:
        # constructor        
        def __init__(self, path):
                self.filePath = path    # path to input file
                self.sequence = list()  # list of integer from filePath
                with open(self.filePath) as f:
                        stringList=f.read().split()
                
                for i in stringList:
                        self.sequence.append(int(i))

        # returns possible maximum sum from a given list 
        def subseqMaxSum(self, sequence):
                initialMax=max(sequence)
                finalMaxSum=currentMaxSum=0

                for i in sequence:
                        currentMaxSum = max(currentMaxSum + i, 0)
                        finalMaxSum = max(finalMaxSum, currentMaxSum)
            
                if initialMax<0 :       # if all elements are less than 0
                        return initialMax
                else:
                        return finalMaxSum

        # returns maximum possible sum of a fixed length (e.g. 5, 10, 20 etc.) out of a longer list
        def subsequenceSumFixedLength(self,n):
                # sanity check
                if len(self.sequence)<n:
                        print("InputError: Third argument should be less than or equal to data size")
                        return -1
                if len(self.sequence)==n:
                        return self.subseqMaxSum(self.sequence)
                if n==1:
                        return max(self.sequence)
                if n<=0:
                        print("InputError: Third argument should be greater than 0")
                        return -1
        
                # computing maximum possible sum within a sublist of n consecutive numbers
                maxSum=tempSum=sum(self.sequence[0:n])
                insideMaxSum=self.subseqMaxSum(self.sequence[0:n])
                maxSum=max(maxSum,insideMaxSum)
                for i in range(1, len(self.sequence)-n+1):
                        tempSum=tempSum+self.sequence[n+i-1]-self.sequence[i-1]
                        insideMaxSum=self.subseqMaxSum(self.sequence[i:n+i])
                        maxSum=max(maxSum,insideMaxSum,tempSum)
                                
                return maxSum

        # computing maximum sum of differences from a list of n integers
        def differenceSum(self,n):
                oriSequence = self.sequence
                sequence = list()

                for i in range(len(oriSequence)-1): #   difference list of each consecutive pairs
                        sequence.append(abs(oriSequence[i]-oriSequence[i+1]))

                if len(oriSequence)<n:  # sanity check
                        print("InputError: Third argument should be less than data size")
                        return -1
                if n<=1:        # sanity check
                        print("InputError: Third argument should be greater than 1")
                        return -1

                n=n-1
                fixedLengthSum=tempSum=sum(sequence[0:n])
                for i in range(1, len(sequence)-n+1):
                        tempSum=tempSum+sequence[n+i-1]-sequence[i-1]
                        fixedLengthSum=max(tempSum,fixedLengthSum)

                return fixedLengthSum


#### Main program
if __name__ == "__main__":
    cwd=os.getcwd()
    path = param = list()
    length = None

    # sanity check and reading input arguments
    if len(sys.argv)==1:
            print("InputError: Path for input data is missing.")
            sys.exit()
    if len(sys.argv)>1:
            path = cwd + "/" + sys.argv[1]
            if os.path.isfile(path)==0:
                    print ("InputError: Input file is not found")
                    sys.exit()
    if  len(sys.argv)>2:
            try:
                    int(sys.argv[2])
            except ValueError:
                    print ("InputError: Third argument should be an integer")
                    sys.exit()
            length = int(sys.argv[2])
    if len(sys.argv)==3:
            print("InputError: Fourth argument is missing. 'values' or 'differences'? ")
    if len(sys.argv)>3:
            param = sys.argv[3]
    if len(sys.argv)>4:
            print("InputError: Arguments should be 2 or 4")

    # executing main tasks    
    if len(sys.argv)==2:
            dataList=SubsequenceMaxSum(path)
            print(dataList.subseqMaxSum(dataList.sequence))
    if len(sys.argv)==4:
            if param=="values":
                    dataList=SubsequenceMaxSum(path)
                    print(dataList.subsequenceSumFixedLength(length))
            elif param=="differences":
                    dataList=SubsequenceMaxSum(path)
                    print(dataList.differenceSum(length))
            else:
                    print("InputError: Fourth argument should be 'values' or 'differences'")


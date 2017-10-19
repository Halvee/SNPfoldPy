import os
import sys
import re
import numpy
import operator
import string
import math
import random
from SequenceProcessing import *

###############################################################################
# This file is part of SNPfold.
#
# SNPfold is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SNPfold is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SNPfold.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


def writefasta(seqname,seq,file_name):#function that writes a seq, seqname and file_name into a fasta file
    ''' writes seqname and seq into a file of name file_name '''
    file = open(file_name,"w")
    file.write('>'+seqname+"\n")
    file.write(seq+'\n')
    file.close()
    return

def process_free_energy_line(line):
    ''' pull out ensemble value from RNAfold output line '''
    line=line.replace(" free energy of ensemble = ","")
    line=line.replace(" kcal/mol","")
    return line

def get_all_ccs(seq1matrix,seq2matrix, minsize=40,cctype="bpProbs",shannonTimeSaver=True):
    ''' extract correlation coefficients for all equivalent comparisons
    of submatrices of size [minsize,minsize]
    to [len(matrix),len(matrix)] '''
    ccArray = numpy.zeros([len(seq1matrix)-minsize,len(seq1matrix)-minsize])

    minCCcurrentsize=minsize
    minCC=1
    minCCi=0

    currentsize=minsize
    while currentsize <= (len(seq1matrix)-minsize):
        i=0
        while(i<len(seq1matrix)-minsize):
            submatrix1 = seq1matrix[i:i+currentsize,i:i+currentsize]
            submatrix2 = seq2matrix[i:i+currentsize,i:i+currentsize]
            if cctype=="shannon":
                cc = numpy.corrcoef(get_shannon(submatrix1,TimeSaver=shannonTimeSaver),get_shannon(submatrix2,TimeSaver=shannonTimeSaver))[0][1]
            else:
                cc = numpy.corrcoef(numpy.sum(submatrix1,axis=0),numpy.sum(submatrix2,axis=0))[0][1]
            if cc < minCC:
                minCC=cc
                minCCi=i
                minCCcurrentsize=currentsize
            ccArray[i,currentsize-minsize]=cc
            i+=1
        currentsize+=1
    return(ccArray,minCC,minCCi,minCCcurrentsize)



def subset_vector(refVector,altVector, positionMatrix):
    filteredRef = []
    filteredAlt = []
    for colIndex in xrange(0,positionMatrix.shape[1]):
        idx = positionMatrix[0,colIndex]
        if refVector[idx]>=0:
            filteredRef.append(refVector[idx])
        if altVector[idx]>=0:
            filteredAlt.append(altVector[idx])
    return refVector,altVector


class RNAsequenceStructureSet:
    ''' takes a list of RNAsequenceStructure instances as input,
    makes it easier to do certain things like get all corrrelation
    coefficients for all variant pairs'''
    def __init__(self,ArrayOfRNAsequenceStructures,constantSeq=True,alignment=None):
        self.cchash={"bpProbs":{'WT':1},"shannon":{"WT":1}}
        self.ccTables={"bpProbs":[["WT",1]],"shannon":[["WT",1]]}
        self.ensembleEhash=dict()
        self.AllFoldingsData = ArrayOfRNAsequenceStructures
        self.globalcc_matrix={"bpProbs":None,"shannon":None}
        self.globalcc_vector={"bpProbs":None,"shannon":None}
        if constantSeq==True:
            self.wt_seq=self.AllFoldingsData[0].SequenceObj.sequence
        self.alignment=None





    def write_as_SNARNASM(self,outfiledir,dataType="BPPROB",seqname="SEQUENCE"):
        if dataType in ["BPPROB","SHANNON"]:
            file = open(outfiledir+dataType+"_data.snarnasm","w")
            file.write("#"+seqname+"\n")
            file.write("#"+self.wt_seq+"\n")
            headers=["SEQPOS"]
            data=[[str(item) for item in xrange(1,len(self.wt_seq)+1)]]
            for item in self.AllFoldingsData:
                headers.append(item.SequenceObj.mutName)
                if dataType=="BPPROB":
                    data.append(list(item.bpProbs))
                elif dataType=="SHANNON":
                    data.append(list(item.shannon))
            data = zip(*data)
            file.write(string.join(headers,"\t")+"\n")
            file.close()
            for row in data:
                file=open(outfiledir+dataType+"_data.snarnasm","a")
                formattedRow = [str(val) for val in row]
                file.write(string.join(formattedRow,"\t")+"\n")
                file.close()
        return self

    def get_ccs(self,dataType="bpProbs"):
        ''' get corrcoeffs between all sequence variants inputted '''
        if type(dataType)==str:
            dataType=[dataType]
        allNames = [Obj.SequenceObj.mutName for Obj in self.AllFoldingsData][1:]
        for mytype in dataType:
            allProbs = [vars(Obj)[mytype] for Obj in self.AllFoldingsData]
            self.globalcc_matrix[mytype]=get_cc(allProbs)
        for mytype in dataType:
            count = 1
            for name in allNames:
                self.cchash[mytype][name]=self.globalcc_matrix[mytype][0,count]
                count+=1
        return self

    def get_pairwise_ccs(self,dataType="bpProbs"):
        ''' get corrcoeffs between all sequence variants inputted '''
        if type(dataType)==str:
            dataType=[dataType]
        for mytype in dataType:
            WTitem = self.AllFoldingsData[0].SequenceObj
            WT = vars(self.AllFoldingsData[0])[myType]
            for Obj in self.AllFoldingsData[1:]:
                MUT = [vars(Obj)[mytype] for Obj in self.AllFoldingsData]
                WT,MUT = subset_vectors(WT,MUT,WTitem.maskedPos)
                allProbs.append(WT)
                allProbs.append(MUT)
            self.globalcc_matrix[mytype]=get_cc(allProbs)
        return self


    def save_metric_set_to_output(self,filename):
        ''' write CC(basepairingProb) and CC(Shannon) to output '''
        header="MUT\tCC_BPPROB\tCC_SHANNON\n"
        file = open(filename,'w')
        file.write(header)
        file.close()
        #print self.cchash
        for SNPobj in self.AllFoldingsData:
            file = open(filename,'a')
            mutName = SNPobj.SequenceObj.mutName

            if mutName!="WT":
                outputLine=[mutName,str(SNPobj.globalcc['bpProbs']['WT']),str(SNPobj.globalcc['shannon']['WT'])]
                file.write("\t".join(outputLine)+"\n")
            file.close()


    def get_P_value_for_ccs(self,ccsOfInterest,subtype="bpProbs"):
        ''' get nonParam. pvalue for a list of CCs '''
        Pvals_for_ccs=[]
        for SNPobj in self.AllFoldingsData:
            row=[SNPobj.mutName,SNPobj.globalcc[subtype]]
            MutTable.append(row)
        sorted(MutTable, key=operator.itemgetter(col))
        for cc in ccsOfInterest:
            index=1
            IsLessThanOrEqualTo = True
            while(index<=len(MutTable) and (IsLessThanOrEqualTo == True)):
                row = MutTable[index]
                if row[1] <= cc:
                    index += 1
                else:
                    IsLessThanOrEqualTo = False
            Pval = float(index)/len(MutTable)
            Pvals_for_ccs.append(Pval)
        return self,Pvals_for_ccs

    def build_mut_table(self,subtype="bpProbs",col=1):
        ''' build and store an (optionally ordered) list of mutations and their corr. coeffs'''
        MutTable=[]
        for SNPobj in self.AllFoldingsData:
            row=[SNPobj.SequenceObj.mutName,SNPobj.globalcc[subtype]['WT']]
            MutTable.append(row)
        if col!=-1:
            MutTable = sorted(MutTable, key=operator.itemgetter(col))
        self.ccTables[subtype]=MutTable
        return self

    def get_P_value_exact(self,SNPsOfInterest,subtype="bpProbs"):
        ''' get exact p-value for SNP of interest '''
        PvalHash=dict()
        nonstandardMuts = []
        MutTable=[]

        for SNPobj in self.AllFoldingsData:
            row=[SNPobj.SequenceObj.mutName,SNPobj.globalcc[subtype]['WT']]
            MutTable.append(row)
        col = 1
        MutTable = sorted(MutTable, key=operator.itemgetter(col))
        self.ccTables[subtype]=MutTable
        rowCount=1
        while (rowCount<=len(MutTable)):
            row = MutTable[rowCount-1]
            if (row[0] in SNPsOfInterest):
                Pval = round(float(rowCount)/(len(MutTable)),5)
                PvalHash[row[0]] = Pval
            rowCount+=1
        for SNP in SNPsOfInterest:
            try:
                PvalHash[SNP]
            except:
                nonstandardMuts.append(SNP)
        return self,PvalHash,nonstandardMuts



class RNAsequenceStructure:
    ''' class that stores data regarding the structure of an instance
    of class SequenceObj from SequenceProcessing.py '''
    def __init__(self,SequenceObj,foldingMethod=None,partitMatrix=None,
             mfeStruct=None,ensembleStruct=None):
        self.SequenceObj = SequenceObj
        self.partitMatrix=partitMatrix
        self.mfeStruct=mfeStruct
        self.mfeE=None
        self.ensembleStruct=None
        self.ensembleE=None
        self.centroidStruct=None
        self.centroidE=None
        self.sampledBinaryStructures=[]
        self.bpProbs=None
        self.shannon=None
        self.bpProbsLogit=None
        self.globalcc={"bpProbs":{"WT":None},"shannon":{"WT":None}}

    def clear_partit_matrix(self):
        del(self.partitMatrix)
        return self

    def get_logit_BpProbs(self):
        ''' calculates and stores log(p/(1-p)),
        given that p = probability of nt. i being basepaired '''
        if self.bpProbs==None:
            return
        self.bpProbsLogit=numpy.zeros(len(self.bpProbs))
        count=0
        while count<len(self.bpProbs):
            self.bpProbsLogit[count]=(math.log((self.bpProbs[count])/(1-self.bpProbs[count])))
            count+=1
        return self


    def fold_RNA(self,typeOfCalc="partit",prog="RNAfold",resultsFull=True,version="2.1.1",additionalOpts="",subsetStart=None,subsetEnd=None):
        ''' folds RNA, can alter commandline opts,
        can fold sequence subset, etc '''
        if (self.SequenceObj.mutantSeq==None):
            seq2fold=self.SequenceObj.sequence
        else:
            seq2fold=self.SequenceObj.mutantSeq

        seq2fold = seq2fold.replace("-","")

        if (typeOfCalc=="partit"):
            if (prog=="RNAfold"):
                foldData = get_RNAfold_partit_matrix(seq2fold,resultsFull=resultsFull,version=version,additionalOpts=additionalOpts)
                if resultsFull==True:
                    (self.partitMatrix,self.mfeStruct,self.mfeE,self.ensembleStruct,self.ensembleE,self.centroidStruct,self.centroidE) = foldData
                    self.bpProbs=numpy.sum(self.partitMatrix,axis=0)
                    if (subsetStart!=None and subsetEnd!=None):
                        self.bpProbs=self.bpProbs[subsetStart:subsetEnd]

        ####OPTIONS FOR CALLING OTHER PROGRAMS CAN GO HERE
        return self

    def get_BpProbsFromMatrix(self,intervalStart=None,intervalEnd=None):
        ''' From partitMatrix, get base pairing probabilities per nucleotide '''
        self.bpProbs = getBpProbPerNt(self.partitMatrix,intervalStart=intervalStart,intervalEnd=intervalEnd)
        return self

    def get_ShannonFromMatrix(self,intervalStart=None,intervalEnd=None,timeSaver=False):
        ''' From partitMatrix, get normalized shannon entropy per nucleotide '''
        self.shannon = get_shannon(self.partitMatrix,intervalStart=intervalStart,intervalEnd=intervalEnd,timeSaver=timeSaver)
        return self




def get_cc(vals,type="perNt",intervalStart=0, intervalEnd=None):
    ''' get CorrelationCoefficient matrix for a set of inputed
    1D vectors (here, usually basepairing probabilities or
    normalized shannon entropies) '''
    count=0
    while count<len(vals):
        #TRIM SEQUENCE IF GETTING CC FOR A SUBSEQUENCES OF INPUTTED VALS
        if intervalEnd==None:
            intervalEnd=len(vals[count])
        vals[count]=vals[count][intervalStart:intervalEnd]
        count+=1
    #vals has to be an array of numpy arrays
    if type=="perNt":
        #vals is assumed to be an array of 1D numpy Arrays
        ccs = numpy.corrcoef(vals)
    elif type=="pairwise":
        count = 1
        ccs = []
        while count<len(vals):
            ccs.append(numpy.corrcoef([vals[count],vals[0]])[0][1])
            count+=1
    return ccs


def getBpProbPerNt(bpProbMatrix,intervalStart=None,intervalEnd=None):
    ''' Get 1D basepairing probabilities vector, given partitMatrix '''
    if intervalStart == None:
        intervalStart=0
    if intervalEnd == None:
        intervalEnd = len(bpProbMatrix)
    # TRIM MATRIX DOWN TO SUBMATRIX IF INTERVALSTART AND/OR INTERVALEND INPUTTED
    MatrixForProbs = bpProbMatrix[intervalStart:intervalEnd,intervalStart:intervalEnd]
    BpProbs=numpy.sum(MatrixForProbs,axis=0)
    return BpProbs

def get_shannon(bpProbMatrix,intervalStart=None,intervalEnd=None,timeSaver=False):
    ''' Get 1D vector of Shannon entropies (normalized against
    the maximum shannon entropy for the same sequence of probabilites, given
    an entirely uniform distribution) given a partitMatrix.   Code adapted from
    Kevin Curran '''
    if intervalStart == None:
        intervalStart=0
    if intervalEnd == None:
        intervalEnd = len(bpProbMatrix)
    MatrixForProbs = bpProbMatrix[intervalStart:intervalEnd,intervalStart:intervalEnd]
    count=0
    nucCount = intervalEnd - intervalStart
    shanValues = numpy.zeros([nucCount])
    #max_se = 0.5*len(bpProbMatrix)*math.log(len(bpProbMatrix),2)
    max_se = math.log(len(bpProbMatrix),2)
    #USE TO SAVE LARGE AMOUNTS OF TIME WHEN GETTING SHANNON ENTROPIES, SLIGHTLY LESS ACCURATE
    if timeSaver==True:
        #approximates zero values to 1E-10, which means shannon values are approximations
        filteredBpProbs = bpProbMatrix.clip(min=1E-10)
        shanValues = -(numpy.sum(filteredBpProbs*numpy.log2(filteredBpProbs),axis=0))/numpy.log2(len(filteredBpProbs))
    else:
        while(count < nucCount):
            shanValue=0
            for p in MatrixForProbs[count]:
                if p > 0:
                    shanValue = shanValue + p*math.log(p,2)
            shanValues[count]= abs(-(1.0/max_se) * shanValue)
            count += 1
    #normalizing factor taken from Kevin Curran's python code
    return shanValues

def get_RNAfold_ensemble_energy(sequence):
    ''' process RNAfold output string to get ensembleEnergy '''
    output = os.popen("echo '"+sequence+"' | RNAfold -p0")
    seqOut = output.readline().rstrip()
    mfe = output.readline().rstrip()
    (mfeStruct,mfeEnergy)=mfe.split(" (")
    mfeEnergy = mfeEnergy.replace(")","")
    mfeEnergy = mfeEnergy.replace("(","")
    ensembleEnergy = output.readline().rstrip()
    mfeFreqInEnsemble = output.readline().rstrip()
    return ensembleEnergy

def get_RNAfold_mfe(sequence):
    ''' process RNAfold output string to get
    mfeStructure, mfe Energy '''
    output = os.popen("echo '"+sequence+"' | RNAfold -p0")
    seqOut = output.readline().rstrip()
    mfe = output.raedline().rstrip()
    (mfeStruct,mfeEnergy)=mfe.split(" ")
    mfeEnergy = mfeEnergy.replace(")","")
    mfeEnergy = mfeEnergy.replace("(","")
    return  mfeStruct,mfeEnergy

def get_RNAfold_partit_matrix(sequence,resultsFull=False,version="2.1.1",additionalOpts="",seqName=None):
    ''' given an input sequence, stores the partitMatrix of basepairing probabilities btwn all
    nucleotides i and j to an instance of numpy.array '''
    nucCount = len(sequence)
    # CREATE RANDOMLY GENERATED OUTFILENAME FOR RNAFOLD WHEN IT IS CALLED
    randFileName=str('.sq'+''.join(random.sample(string.letters+string.digits,6)))
    file_name=randFileName+"_dp.ps"
    extra_file_name = randFileName+"_ss.ps"
    # WRITE SEQ TO OUTFILE
    writefasta(randFileName,sequence,randFileName+'.fa')
    # CALL RNAFOLD.. startSearchTerm IS BEGINNING OF WHERE VALUES ARE LISTED IN POSTSCRIPT OUTFILE
    output = os.popen('RNAfold -p '+additionalOpts+' <'+randFileName+'.fa')#runs RNAfold
    #os.system("rm "+extra_file_name)
    #PROCESS RNAFOLD STDOUT DATA
    header = output.readline().rstrip()
    seqOut = output.readline().rstrip()
    mfe = output.readline().rstrip()
    mfeStruct = mfe[0:nucCount]
    mfeE = mfe[nucCount+1:]
    mfeE = mfeE.replace("(","")
    mfeE = mfeE.replace(")","")
    mfeE = float(mfeE)
    ensemble = output.readline().rstrip()
    ensembleStruct = ensemble[0:nucCount]
    ensembleE = ensemble[nucCount+1:]
    ensembleE = ensembleE.replace("[","")
    ensembleE = ensembleE.replace("]","")
    ensembleE = float(ensembleE)
    centroid = output.readline().rstrip()
    centroidStruct = centroid[0:nucCount]
    centroidE = centroid[nucCount+1:]
    centroidE = centroidE.replace("{","")
    centroidE = centroidE.replace("}","")

    # INITIALIZE A NUMPY ARRAY ZEROS OF DIM [nucCount,nucCount]
    myMatrix = numpy.zeros([nucCount,nucCount])

    # READ IN POSTSCRIPT OUTPUT CONTAINING PARTITMATRIX
    file = open(file_name,"r")
    results_data = file.read()
    file.close()

    os.system("rm "+randFileName+"*")
    #FIND WHERE THE PERTINENT DATA IS IN THE POSTSCRIPT FILE
    try:
        startSearchTerm = '%start of base pair probability data'
        string.index(results_data,startSearchTerm)
    except:
        startSearchTerm = '%data starts here'
        string.index(results_data,startSearchTerm)
    startLocation=string.find(results_data,startSearchTerm)
    stopLocation=string.find(results_data,'showpage')
    results_data = results_data[(startLocation+len(startSearchTerm)+1):stopLocation]
    rows = results_data.splitlines()

    for row in rows:
        #FOR EACH ROW, ID THE X POS, Y POS, SQRT BPPROB VALUE
        values=row.split()
        x=values[0]
        y=values[1]
        partit=values[2]
        # If the line had 'ubox' in it (which means that it is a part of the partit. func.):
        # Fill in at [position x, position y] of the matrix the value partit^2
        # Fill in at [position y, position x] of the matrix the value partit^2
        if values[3]=='ubox':
            myMatrix[int(y)-1,int(x)-1]=(float(partit)*float(partit))
            myMatrix[int(x)-1,int(y)-1]=(float(partit)*float(partit))


    # These values should be the same, since the matrix is Hermitian
    # Need s for sum of values in columns
    s=numpy.sum(myMatrix,axis=1)
    s=numpy.sum(myMatrix,axis=0)

    # Store sums in the list named 'partit'
    if resultsFull == True and seqName == None:
        return myMatrix,mfeStruct,mfeE,ensembleStruct,ensembleE,centroidStruct,centroidE
    elif resultsFull == True and seqName !=None:
        return myMatrix,mfeStruct,mfeE,ensembleStruct,ensembleE,centroidStruct,centroidE,seqName
    elif resultsFull == False and seqName!=None:
        return myMatrix,seqName
    else:
        return myMatrix

def get_RNAfold_partit_matrix_MultiThreaded(sequence,variantName,additionalOpts):
    ''' A function involved in running RNAfold in a multithreaded fashion.   Functions that carry out the
    multithreading are in SNPfold_commandline.py'''
    if variantName==None:
        variantName="WT"
    nucCount = len(sequence)
    randFileName=str('.sq'+''.join(random.sample(string.letters+string.digits,6)))
    file_name=randFileName+"_dp.ps"
    writefasta(randFileName,sequence,randFileName+'.fa')
    output = os.popen('RNAfold -p '+additionalOpts+' < '+randFileName+'.fa')

    #calculate the matrix of base-pairing probabilities for a given sequence
    header = output.readline().rstrip()
    seqOut = output.readline().rstrip()
    mfe = output.readline().rstrip()
    mfeStruct = mfe[0:nucCount]
    mfeE = mfe[nucCount+1:]
    mfeE = mfeE.replace("(","")
    mfeE = mfeE.replace(")","")
    mfeE = float(mfeE)
    ensemble = output.readline().rstrip()
    ensembleStruct = ensemble[0:nucCount]
    ensembleE = ensemble[nucCount+1:]
    ensembleE = ensembleE.replace("[","")
    ensembleE = ensembleE.replace("]","")
    ensembleE = float(ensembleE)
    centroid = output.readline().rstrip()
    centroidStruct = centroid[0:nucCount]
    centroidE = centroid[nucCount+1:]
    centroidE = centroidE.replace("{","")
    centroidE = centroidE.replace("}","")

    #get nucCount from the length of the RNA sequence
    #create a matrix of zeros of RNA length X RNA length
    myMatrix = numpy.zeros([nucCount,nucCount])

    #get partition function matrix data from RNAfold output
    file = open(file_name,"r")
    results_data = file.read()
    file.close()
    try:
        startSearchTerm = '%start of base pair probability data'
        string.index(results_data,startSearchTerm)
    except:
        startSearchTerm = '%data starts here'
        string.index(results_data,startSearchTerm)


    os.system("rm "+randFileName+"*")

    #format the results data from the output file... find where the pertinent data is
    startLocation=string.find(results_data,startSearchTerm)
    stopLocation=string.find(results_data,'showpage')
    results_data = results_data[(startLocation+len(startSearchTerm)+1):stopLocation]
    rows = results_data.splitlines()
    for row in rows:
        #for every row, identify the x position, y position and partit value
        values=row.split()
        x=values[0]
        y=values[1]
        partit=values[2]
        # If the line had 'ubox' in it (which means that it is a part of the partit. func.:
        # Fill in at [position x, position y] of the matrix the value partit^2
        # Fill in at [position y, position x] of the matrix the value partit^2
        if values[3]=='ubox':
            myMatrix[int(y)-1,int(x)-1]=(float(partit)*float(partit))
            myMatrix[int(x)-1,int(y)-1]=(float(partit)*float(partit))


    # These values should be the same, since the matrix is hermitian
    # Need s for sum of values in columns
    s=numpy.sum(myMatrix,axis=1)
    s=numpy.sum(myMatrix,axis=0)

    # Store sums in the list named 'partit'
    return variantName,myMatrix,mfeStruct,mfeE,ensembleStruct,ensembleE,centroidStruct,centroidE




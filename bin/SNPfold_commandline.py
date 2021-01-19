#!/usr/bin/env python

###############################################################################
# This file is part of SNPfoldPy.
#
# SNPfoldPy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SNPfoldPy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SNPfold.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################



############################################################
#
# UPDATE LOG
#
# V1.03 (01/19/2021)
# * Output directory path now written to standard error.
# * fixed bug where runs without accurate p-value turned on weren't having
#   results saved to file even though option '-s' was being passed
# * fixed bug where program would try and fail to run when trying to calculate
#   accurate p-values for multi-nucleotide mutants. Accurate p-value algorithm
#   described in original SNPfold paper is only valid for a single nucleotide
#   variant within the input sequence, and not for insertions, deletions, or
#   alterations of more than one nucleotide at a time.
# * full help descriptions printed in scenario where positional arguments not
#   fully provided
# * slight expansion to option help descriptions
# V1.02 (10/18/2017)
# * Changed formal name of package to SNPfoldPy and put it up on github.
#   Changed all references of package in documentation to SNPfoldPy.
# * Fixed memory overlap issue by
# V1.01 (6/26/2013)
# * In 'detect indels' function from SequenceProcessing.py, changed regular
#   expression search to search all subSNPs for a given mutant for positions
#   in center (ie. 'A44G,G50C')
# * Made the same change to 'detect_strechPtMutations' function in
#   SequenceProcessing.py
# * Added a 'mutantSeq' variable initialization in the
#   __init__ for class SequenceObj in SequenceProcessing.py
# * In 'read_data_mutants', set up a while loop to remove empty lines from
#   mutation input files
#
############################################################




import os
import sys
import re
import math
import numpy
import random
import time
from shutil import rmtree
from multiprocessing import Pool,freeze_support


from SNPfold.BasicFunctions import *


''' global variables '''
import time
MultiThreadingJobs=[]
percentDone = 0


def mainpulate_outfile(filename,append2=None):
    ''' create empty file, or append to pre-existing '''
    if append2==None:
        file=open(filename, "w")
        file.write("")
        file.close()
    else:
        file =open(filename, "a")
        file.write(append2+"\n")
        file.close()
    return

def update_percent_done(currentCount,totalCount,percentDoneByIncrem,increment=0.1):
    ''' update the percent completion of the task being run '''
    global percentDone
    currentCount+=1
    ratio = float(currentCount)/totalCount
    percentDone = ratio
    if percentDoneByIncrem+increment <= percentDone:
        percentDoneByIncrem+=increment
        myPercentDone = percentDoneByIncrem * 100
        sys.stderr.write(str(myPercentDone)+"% completed\n")
    return percentDoneByIncrem

def multithreading_log_result(result):
    ''' store result variable to MultiThreadingJobs global variable '''
    global MultiThreadingJobs
    MultiThreadingJobs.append(result)

def generate_random_string(stringLength):
    ''' generate a string of eight random characters,
    pulled from a string of possible characters '''
    allChars=(string.ascii_letters+string.digits)
    randomString=''
    count=0
    while len(randomString)<stringLength:
        newChar=random.choice(allChars)
        randomString+=newChar
    return randomString


def read_data(filename,asFasta=False):
    ''' read in data from a specified file '''
    file=open(filename,'r')
    file_contents=file.read()
    file.close()
    if (asFasta==True):
        try:
            fileText=file_contents.splitlines()
            seqInfo=fileText[0]
            seqName = seqInfo.replace(">","")
            sequence="".join(fileText[1:])
            return seqName,sequence
        except:
            usage()
            sys.exit(2)
    return file_contents

def read_data_mutants(filename):
    file=open(filename,'r')
    file_contents=file.readlines()
    file.close()
    file_contents = [line.rstrip() for line in file_contents]
    while '' in file_contents:
        file_contents.remove('')
    return ":".join(file_contents)

def direxists(dirName):
    ''' gives user option to overwrite output directory if it already exists '''
    print 'That directory already exists. Do you want to overwrite it?\n\
Everything will be deleted. Type "y" to overwrite, "n" to cancel.'
    overwrite = raw_input("> ")

    if overwrite in ['y','Y']:
        rmtree(os.getcwd()+'/output/'+dirName)
        print '%s was emptied.'%(os.getcwd()+'/output/'+dirName)
        os.system('mkdir output/'+dirName)
        sys.stderr.write('output dir : %s' %(os.getcwd()+'/output/'+dirName+"\n"))
    elif overwrite in ['n','N']:
        sys.exit('Operation aborted!\n')
    else:
        print 'Invalid input.'
        direxists(dirName)

def SequenceSanityCheck(seq):
    ''' checks to see if nucleotides in sequence are normal RNA nucleotides '''
    ntsAllowed = ['A','C','G','T','U']
    seqIndex = 0
    while(seqIndex<len(self.sequence)):
        if (sequence[seqIndex] not in self.ntsAllowed):
            print "character ("+str(seq[seqIndex])+") not allowed in folding"
            sys.exit()
            return False
        seqIndex+=1
    return True


def check_for_errors(wild_type, SNPs):
    ''' filter for mutation formatting errors '''
    mutants = []
    bad_SNPs = []
    nts = ['A','T','U','G','C']

    SequenceSanityCheck(wild_type)
    for SNP in SNPs:
        errors=False
        subSNPs=SNP.split(",")
        for subSNP in subSNPs:
            positionset=re.findall(r"\d+",subSNP)
            if len(positionset)==1:
                position=positionset[0]
                pos=int(position)-1
                [wt,mut]=string.split(subSNP,position)
                if SequenceSanityCheck(wt) == False:    errors = True
                if SequenceSanityCheck(mut) == False:   errors = True
                if (wt_nt != wild_type[pos:pos+len(wt_nt)]):    errors=True
            else:
                errors=True

        if errors==False:
            mutants.append(SNP)
        else:
            bad_SNPs.append(SNP)


    if len(bad_SNPs)>0:
        print "ERROR WITH THE FOLLOWING SNPs : "
        for item in bad_SNPs:
            print item
        sys.exit(2)

    return

def process_args(args):
    ''' handles conflicting user input args '''
    if os.path.isfile(args.wild_type)==True:
        args.seqName,args.wild_type=read_data(args.wild_type,asFasta=True)
    if os.path.isfile(args.mutants)==True:
        try:
            args.mutants=read_data_mutants(args.mutants)
        except:
            sys.exit(2)

    elif args.mutants=="ALL":
        wtSeqObj = SequenceObj(args.wild_type)
        args.mutants=string.join(wtSeqObj.allPossiblePointMuts(),":")



    if os.path.isdir('output')==False:
        os.system('mkdir output')
    if  detect_indels(args.mutants.split(":"))==True:
        print "error : no insertion or deletion mutations"
        sys.exit(1)
    if detect_strechPtMutations(args.mutants.split(":"))==True:
        args.accurate=False

    if args.accurate==True or args.nameDir!=None or args.SNARNASM==True or args.save==True:
        args.save=True
        if args.nameDir==None:
            args.nameDir=generate_random_string(8)
        if os.path.isdir('output/'+args.nameDir)==True:
            direxists(args.nameDir)
        else:
            os.system('mkdir output/'+args.nameDir)
            sys.stderr.write('output dir : %s' %(os.getcwd()+'/output/'+args.nameDir+"\n"))

    if args.additionalOpts!="":
        args.additionalOpts = "-"+args.additionalOpts
        if args.RNAfoldVersion=="2.1.1":
            args.additionalOpts = "-"+args.additionalOpts


    args.wild_type=str.upper(args.wild_type)



    args.mutants=str.upper(args.mutants)
    args.mutants = DNAtoRNA(args.mutants)

    return args

def SNPfold_commandline(argv):
    ''' runs SNPfold from commandline '''
    import argparse
    #check to see if output dir exists.. if it doesn't, it is created
    #establish acceptable nts in seq and mutations listed


    percentDoneByIncrem = 0.0

    opts = argparse.ArgumentParser()
    opts.add_argument("-A","--accurate",
                      help="Calculate accurate p-values for mutations " + \
                           "indicated by user. Algorithm only works with " + \
                           "single nucleotide mutants.",
                      action="store_true")
    opts.add_argument("-m","--metric",choices=['bpProbs','cc','shannon','all'],
                      default='cc',
                      help="Specify metrics to be written to file, in the " + \
                           "scenario that the '--save' option is turned on")
    opts.add_argument("-s","--save",action="store_true",
                      help="results are saved to file(s) rather than just " + \
                           "being printed to standard output.")
    opts.add_argument("-n","--nameDir",
                      help="Name of the directory where output files " + \
                           "would be written to.")
    opts.add_argument("-t","--threads",
                      help="if calculating accurate p-values, number of processors to use for threading",
                      type=int,default=1)
    opts.add_argument('-a','--additionalOpts',
                      help="additional opts for RNAfold to use",
                      type=str,default="")
    opts.add_argument('--keep-partit-matrices',
                      help="keep full matrices in memory instead of deleting",
                      action="store_true", default=False)
    opts.add_argument("-S","--SNARNASM",
                      help="output data in SNARNASM format",
                      action="store_true")
    opts.add_argument("-N","--seqName",help="name of sequence",
                      default="SEQUENCE")
    opts.add_argument("wild_type", 
                      help="Nucleotide sequence. Should only " + \
                           "consist of the following nucleotides : " + \
                           "A, C, G, U, T."
                     )
    opts.add_argument("mutants", 
                      help="SNPs separated by colon, e.g. A5G:T18C . " + \
                           "Multi-nucleotide mutants can be passed " + \
                           "in a comma-seperated set of single nucleotide " + \
                           "mutations (ex: C1G,G3U:A5G:T18C) " + \
                           "but the accurate p-value algorithm cannot be " + \
                           "used on these mutants since it is meant " + \
                           "for single nucleotide mutants only.")
    if len(argv) < 3: 
        opts.print_help()
        print("\nerror : too few arguments")
        print("SNPfold_commandline.py [OPTS] <wild_type> <mutants>\n")
        sys.exit(1)
    args = opts.parse_args()
    args = process_args(args)

    wtSeqObj = SequenceObj(args.wild_type,mutation=None,
                           seqName=args.seqName,mutName="WT")

    allVariants=[wtSeqObj]
    AllSNPsToParse=None
    SNPs=args.mutants.split(":")
    #check_for_errors(args.wild_type,SNPs)
    if args.accurate == False:
        MutSNPsToParse = args.mutants.split(":")
        SNPsToParse = MutSNPsToParse
    else:
        if args.mutants.find(",") != -1:
            sys.exit("error : multi-nucleotide mutants not " + \
                     "compatible with accurate p-value calculation.")
        MutSNPsToParse = args.mutants.split(":")
        SNPsToParse = wtSeqObj.allPossiblePointMuts()

    for SNP in SNPsToParse:
        mutSeqObj = SequenceObj(args.wild_type,mutation=SNP,seqName=args.seqName,mutName = SNP)
        allVariants.append(mutSeqObj)
    allVariantsStructures=[]
    count = 0
    for variant in allVariants:
        ''' create and store RNAsequenceStructure instances '''
        StructureObj = RNAsequenceStructure(variant)
        allVariantsStructures.append(StructureObj)

    if args.threads<=1:
        ''' no multiprocessing '''
        count = 0
        ''' fold each RNA, get bpProbs, normalized Shannon entropy '''
        for StructureObj in allVariantsStructures:
            StructureObj.fold_RNA(resultsFull=True,additionalOpts=args.additionalOpts)
            StructureObj.get_ShannonFromMatrix(timeSaver=False)
            if args.keep_partit_matrices == False:
                StructureObj.clear_partit_matrix()
            if args.accurate==True:
                percentDoneByIncrem = update_percent_done(count,len(allVariants),percentDoneByIncrem,increment=0.1)
            count+=1

    else:
        ''' multiprocessing '''
        p = Pool(processes=args.threads)
        allVariantsStructures[0].SequenceObj.mutation="WT"
        p.apply_async(get_RNAfold_partit_matrix_MultiThreaded,
                      args = (allVariantsStructures[0].SequenceObj.sequence,
                              allVariantsStructures[0].SequenceObj.mutation,
                              args.additionalOpts),
                      callback=multithreading_log_result)
        count=1
        while count<len(allVariantsStructures):
            p.apply_async(get_RNAfold_partit_matrix_MultiThreaded,
                          args = (allVariantsStructures[count].SequenceObj.mutantSeq,
                                  allVariantsStructures[count].SequenceObj.mutation,
                                  args.additionalOpts),
                          callback=multithreading_log_result)
            count+=1

        p.close()
        # percent done counter goes here
        p.join()


        global MultiThreadingJobs
        MultiThreadingResults=dict()
        for result in MultiThreadingJobs:
            MultiThreadingResults[result[0]]=result
        count = 0
        ''' store structure predictions to RNAsequenceStructure instances '''
        while count<len(allVariantsStructures):
            variantName=allVariantsStructures[count].SequenceObj.mutation
            allVariantsStructures[count].partitMatrix=MultiThreadingResults[variantName][1]
            allVariantsStructures[count].mfeStruct=MultiThreadingResults[variantName][2]
            allVariantsStructures[count].mfeE=MultiThreadingResults[variantName][3]
            allVariantsStructures[count].ensembleStruct=MultiThreadingResults[variantName][4]
            allVariantsStructures[count].ensembleE=MultiThreadingResults[variantName][5]
            allVariantsStructures[count].centroidStruct=MultiThreadingResults[variantName][6]
            allVariantsStructures[count].centroidE=MultiThreadingResults[variantName][7]
            allVariantsStructures[count].bpProbs=numpy.sum(MultiThreadingResults[variantName][1],axis=0)
            allVariantsStructures[count].shannon=get_shannon(allVariantsStructures[count].partitMatrix)
            if args.keep_partit_matrices == False:
                del(allVariantsStructures[count].partitMatrix)

            count+=1



    AllFoldingsObj = RNAsequenceStructureSet(allVariantsStructures)
    AllFoldingsObj.get_ccs(dataType=["bpProbs","shannon"])
    count=1
    while count<len(allVariantsStructures):
        # collect cc btwn wt, each mut
        allVariantsStructures[count].globalcc["bpProbs"]["WT"]=AllFoldingsObj.globalcc_matrix["bpProbs"][0,count]
        allVariantsStructures[count].globalcc["shannon"]["WT"]=AllFoldingsObj.globalcc_matrix["shannon"][0,count]
        count+=1

    if args.accurate==False:
        if args.metric=='cc' or args.metric=="all":
            ''' default, returns only correlation coefficient btwn wildtype, each mut '''
            res=["\t".join(["MUT","CC_BPPROB","CC_SHANNON"])]

            count = 1

            while count<len(allVariantsStructures):
                # collect cc btwn wt, each mut
                sequenceData = allVariantsStructures[count].SequenceObj
                if sequenceData.mutName!="WT" and args.accurate != True:
                    res.append("\t".join([sequenceData.mutation,
                                          str(allVariantsStructures[count].globalcc["bpProbs"]["WT"]),
                                          str(allVariantsStructures[count].globalcc["shannon"]["WT"])]))
                count+=1

            ''' if desired by user, write to file, otherwise print to stdout '''
            if args.save==True:
                mainpulate_outfile("output/"+args.nameDir+"/input_mutations_results.txt")
                for str_i in res:
                    mainpulate_outfile("output/"+args.nameDir+"/input_mutations_results.txt",
                                       append2=str_i)
            else:
                for str_i in res: print(str_i)


    else:
        ''' calculate, report p-values for each mutation of interest '''
        if args.save==True:
            mainpulate_outfile("output/"+args.nameDir+"/input_mutations_allResults.txt")
            ''' save all possible point mutations correlation coefficients '''
            for subtype in allVariantsStructures[0].globalcc.keys():
                i = 1
                AllFoldingsObj.build_mut_table(subtype=subtype,col=-1)
                mainpulate_outfile("output/"+args.nameDir+"/all_pointmuts_"+subtype+"_CCs.txt")
                mainpulate_outfile("output/"+args.nameDir+"/all_pointmuts_"+subtype+"_CCs.txt",
                                    append2="\t".join(["MUT",string.upper(subtype)]))
                while i < len(AllFoldingsObj.ccTables[subtype]):
                    AllFoldingsObj.ccTables[subtype][i][1]=str(AllFoldingsObj.ccTables[subtype][i][1])
                    if args.save==True:
                        mainpulate_outfile("output/"+args.nameDir+"/all_pointmuts_"+subtype+"_CCs.txt",
                                            append2="\t".join(AllFoldingsObj.ccTables[subtype][i]))
                    i+=1






        PvalSets=dict()
        for subtype in allVariantsStructures[0].globalcc.keys():
            x = AllFoldingsObj.get_P_value_exact(SNPs,subtype=subtype)
            Pvals = x[1]
            PvalSets[subtype]=Pvals
            NonStandardMuts = x[2]

        print "\t".join(["MUT","CC_BPPROB","CC_BPPROB_PVAL","CC_SHANNON","CC_SHANNON_PVAL"])
        mainpulate_outfile("output/"+args.nameDir+"/input_mutations_allResults.txt",
                            append2="\t".join(["MUT","CC_BPPROB","CC_BPPROB_PVAL","CC_SHANNON","CC_SHANNON_PVAL"]))


        #print PvalSets[subtype].keys()

        for mut in SNPs:
            outputLine=[mut]
            for subtype in AllFoldingsObj.globalcc_matrix.keys():
                outputLine.append(str(AllFoldingsObj.cchash[subtype][mut]))
                outputLine.append(str(PvalSets[subtype][mut]))
            print "\t".join(outputLine)
#506
            if args.save==True:
                mainpulate_outfile("output/"+args.nameDir+"/input_mutations_allResults.txt"
                                ,append2="\t".join(outputLine))

    ''' determine output file names based on whether or not accurate pvals used '''
    fprefix="input_mutations"
    if args.accurate == True:
        fprefix="all_pointmuts"

    if args.metric=='bpProbs' or args.metric=="all":
        ''' output bpProbs per sequence variant '''
        bpProbsTable=[Obj.bpProbs for Obj in allVariantsStructures]
        if args.save==True:
            i = 0
            mainpulate_outfile("output/"+args.nameDir+"/"+fprefix+"_bpProbs_perNt.txt")
            while i < len(AllFoldingsObj.AllFoldingsData):
                bpProbs = AllFoldingsObj.AllFoldingsData[i].bpProbs
                if args.save==True:
                    output = [AllFoldingsObj.AllFoldingsData[i].SequenceObj.mutName]
                    output.extend([str(val) for val in bpProbs])
                    output = "\t".join(output)
                    mainpulate_outfile("output/"+args.nameDir+"/"+fprefix+"_bpProbs_perNt.txt",
                                        append2=output)

                i+=1


    if args.metric=='shannon' or args.metric=="all":
        ''' output normalized shannon ent. per sequence variant '''
        shannonTable=[Obj.shannon for Obj in allVariantsStructures]
        if args.save==True:
            i = 0
            mainpulate_outfile("output/"+args.nameDir+"/"+fprefix+"_shannon_perNt.txt")
            while i < len(AllFoldingsObj.AllFoldingsData):
                shannon = AllFoldingsObj.AllFoldingsData[i].shannon
                if args.save==True:
                    output = [AllFoldingsObj.AllFoldingsData[i].SequenceObj.mutName]
                    output.extend([str(val) for val in shannon])
                    output = "\t".join(output)
                    mainpulate_outfile("output/"+args.nameDir+"/"+fprefix+"_shannon_perNt.txt",
                                        append2=output)
                i+=1


    if args.SNARNASM==True:
        AllFoldingsObj.write_as_SNARNASM("output/"+args.nameDir+"/",dataType="BPPROB",seqname=args.seqName)
        AllFoldingsObj.write_as_SNARNASM("output/"+args.nameDir+"/",dataType="SHANNON",seqname=args.seqName)

    return


if (__name__=="__main__"):
    # ADJUSTS INPUT ACCORDINGLY IF INPUT MUT STARTS WITH '-'
    if sys.argv[-1][0]=="-":
        sys.argv.insert(-1,"--")
    SNPfold_commandline(sys.argv)

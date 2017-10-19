import os
import sys
import sys
import re
import string
import tempfile
import numpy

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

def reverseComp(seq,type="RNA"):
	''' gets reverse complement of input sequence '''
	seq = seq.translate(string.maketrans("AGCTUagctu", "TCGAAtcgaa"))
	seq = seq[::-1]
	# converts sequence to RNA, if specified by user
	if type=="RNA":	seq=DNAtoRNA(seq)
	return seq

def DNAtoRNA(seq,changeToCaps=True):
	''' covnerts DNA sequence to RNA sequence '''
	try:
		if changeToCaps==True:
			seq = seq.upper()
		if seq != None:
			seq = seq.replace("T","U")
			seq = seq.replace("t","u")
		
	except:
		pass
	return seq

def merge_intervals(intervalSet):
	intervalSet.sort()
	newSet = []
	count = 1
	while count < len(intervalSet):
		if ((intervalSet[count][0] > intervalSet[count-1][0]) and (intervalSet[count][0] <= intervalSet[count-1][1])):
			intervalSet[count-1] = (intervalSet[count-1][0],max(intervalSet[count-1][1],intervalSet[count][1]))
			del(intervalSet[count])
		else:
			count+=1
	return intervalSet
			   
		
		

def get_interval(chrom,chromstart,chromend,strand="+",species="hg19",path=os.getenv("HOME")+"/genomes/",toRNA=False):
	'''given the possession of fasta files containing chromosome sequences
	from the reference genome in question, extracts the sequence
	at the found at the input genomic coordinates'''
	#NOTE : REQUIRES SAMTOOLS TO BUILD SEQUENCE INDECES
	chromstart=int(chromstart)
	chromend=int(chromend)
	pathToGenome=path+species+"/"
	filename=pathToGenome+"hs_ref_"+chrom+".fa"
	filenameIndex=filename+".fai"
	if os.path.exists(filenameIndex)==False:
		os.system("samtools faidx "+filename)
	#ONLY ONE LINE EXPECTED IN EACH CHROMOSOME INDEX FILE
	fileINindex=open(filenameIndex,"r")
	line = fileINindex.readline()
	fileINindex.close()
	line=line.rstrip()
	[fastaheaderstuff,bytelength,bytestart,linebytelength,endoflinebytepos]=line.split("\t")
	''' store the byte length, start byte position to read from, 
	length of each line in file in bytes, location of line end in bytes '''
	bytelength=int(bytelength)
	bytestart=int(bytestart)
	linebytelength=int(linebytelength)
	endoflinebytepos=int(endoflinebytepos)
	fastafile=open(filename,"r")

	#CALCULATE THE START AND END BYTE POSITION OF INTERVAL
	numberOfNewLinesBefore=int(chromstart)/int(linebytelength)
	numberOfNewLinesBeforeEnd=int(chromend)/int(linebytelength)
	byteStartPoint=bytestart+chromstart+numberOfNewLinesBefore
	byteEndPoint=bytestart+chromend+numberOfNewLinesBeforeEnd
	intervalLength=int(byteEndPoint)-int(byteStartPoint)

	#EXTRACT INTERVAL
	fastafile.seek(int(byteStartPoint))
	sequence=fastafile.read(intervalLength)
	fastafile.close()
	sequence=sequence.replace("\n","")
	if strand=="-":
		sequence=reverseComp(sequence,type="DNA")
	if toRNA==True:
		sequence=DNAtoRNA(sequence)
	return sequence

def changeTolocalCoord(startcoord,endcoord,relToStart,relToEnd,strand):
	''' Given a pair of stop/start coordinates in genome, and a 
	pair of coordinates interior to the first pair of coords,
	find the local start/stop coordinates for the interior pair'''
	if strand=="+":
		startcoord=startcoord-relToStart
		endcoord=endcoord-relToStart
	elif strand=="-":
		_startcoord=-(endcoord-relToEnd)
		_endcoord=-(startcoord-relToEnd)
		startcoord=_startcoord
		endcoord=_endcoord
	return startcoord,endcoord

def SequenceObjFromCoords(chrom,start,stop,name=None,strand=None,type="hg19"):
	''' Get an instance of SequenceObj from interval data '''
	sequence=get_interval(chrom,start,stop,strand=strand)
	MySequenceObj=SequenceObj(sequence,seqName=name,seqChrom=chrom,seqChromStart=start,seqChromEnd=stop)
	return MySequenceObj

def generate_allPossiblePointMuts(sequence):
	''' Generate all possible point mutations for the sequence '''
	mutations = []
	allNts='ACGU'
	position = 0
	while (position < len(self.sequence)):
		nt = self.sequence[position]
		NtChoices = allNts.replace(nt,"")
		for NtChoice in NtChoices:
			mutation = nt+str(position+1)+NtChoice
			mutations.append(mutation)
		position+=1
	return mutations

def sizeNormalizeSeqs(seq1,seq2):
	while(1):
		if len(seq1)>len(seq2):
			seq2 = seq2 + "-"
		elif len(seq2)>len(seq1):
			seq1 = seq1 + "-"
		else:
			break
	return seq1,seq2

def detect_strechPtMutations(SNPs):
	''' detect point mutations more than one nucleotide in length '''
	stretchMuts = False
	for SNP in SNPs:
		subSNPs = SNP.split(",")
		for subSNP in subSNPs:
			positionset=re.findall(r"\d+",subSNP)
			wt,mut = subSNP.split(positionset[0])
			if len(wt) != 1 or len(mut)!=1:
				stretchMuts = True
	return stretchMuts

def detect_indels(SNPs):
	indels = False
	for SNP in SNPs:
		subSNPs = SNP.split(",")
		for subSNP in subSNPs:
			positionset=re.findall(r"\d+",subSNP)
			wt,mut = subSNP.split(positionset[0])
			if "-" in (wt,mut):
				indels = True
	return indels
	
class SequenceObj(object):
	''' Class for a sequence extracted from some genomic interval,
	will allow you to get local coord for mutation in interval, 
	given genomic coords for interval start/stop, mutation start/stop'''
	def __init__(self,sequence,mutation=None,seqName=None,seqChrom=None,seqChromStart=None,\
seqChromEnd=None, mutChromStart=None, mutChromEnd=None, mutName=None, strand="None",\
refNts=None,altNts=None,localCoordStart=None,localCoordEnd=None,sameStrand=False,assumeWTseq=True,\
mutantSeq=None):
		sequence = DNAtoRNA(sequence)
		mutation = DNAtoRNA(mutation)
		refNts = DNAtoRNA(refNts)
		altNts = DNAtoRNA(altNts)
		sequence=string.upper(sequence)
		self.sequence = sequence
		self.mutation = mutation
		self.seqName = seqName
		self.seqChrom = seqChrom
		self.seqChromStart = seqChromStart
		self.seqChromEnd = seqChromEnd
		self.mutName = mutName
		self.strand = strand
		self.mutChromStart = mutChromStart
		self.mutChromEnd = mutChromEnd
		self.refNts=refNts
		self.altNts=altNts
		self.localCoordStart=localCoordStart
		self.localCoordEnd=localCoordEnd
		self.sameStrand=sameStrand
		self.Substitutions=[]
		self.Insertions=[]
		self.Deletions=[]
		self.insOffset=0
		self.mutantSeq=mutantSeq

		self.SNPPositions=[]
		self.InsPositions=[]
		self.DelPositions=[]
		
		self.pairwiseAln = None
		self.pairwisePos = None
		self.maskedAln = None
		self.ntsAllowed = ['A','C','G','T','U']
		self.mutCharsAllowed=['A','C','G','T','U','-']
		self.badSNP=False

		self.SequenceSanityCheck()
		
		if self.seqChromStart!=None and self.seqChromEnd!=None:
			if self.mutChromStart!=None and self.mutChromEnd!=None:
				if self.refNts!=None and self.altNts!=None:
					# CHANGE MUTATION TO LOCAL MUTATION
					self.changeToLocalMutation()
					if self.strand!="None" and sameStrand==False:
						self.changeNtsBasedOnStrand()
					self.localCoordStart,self.localCoordEnd=changeTolocalCoord(self.mutChromStart,self.mutChromEnd,self.seqChromStart,self.seqChromEnd,self.strand)
					# GIVE NAME TO MUTATION ON THE INTERVAL LEVEL
					self.mutation=self.refNts+str(self.localCoordEnd)+self.altNts
		if self.mutation!=None:
			self.ProcessMutation()
		elif self.mutName=="WT":
			self.mutation="WT"

	def mask_aln(self):
		''' mask positions in wt/mut alignment where structural info can't be compared '''
		self.maskedAln = numpy.array([self.pairwiseAln[0],self.pairwiseAln[2]])
		self.maskedPos = numpy.array([self.pairwiseAln[1],self.pairwiseAln[3]])
		colIndex = 0
		#print "WT SEQ : "
		#print string.join(self.pairwiseAln[0],"")
		#print "MUT SEQ : "
		self.mutantSeq=string.join(self.pairwiseAln[2],"")
		for colIndex in xrange(0,self.maskedAln.shape[1]):
			if self.maskedAln[0,colIndex]=="-" and self.maskedPos[0,colIndex]==-1:
				self.maskedAln[1,colIndex]="*"
				self.maskedPos[1,colIndex]=-2
			elif self.maskedAln[1,colIndex]=="-" and self.maskedPos[1,colIndex]==-1:
				self.maskedAln[0,colIndex]="*"
				self.maskedPos[0,colIndex]=-2
		return self

	def changeToLocalMutation(self):
		''' calculate local coords for mutation based on genomic coords '''
		if self.strand=="+":
			self.loocalCoordStart=self.mutChromStart-self.seqChromStart
			self.loocalCoordEnd=self.mutChromEnd-self.seqChromStart
		elif self.strand=="-":
			self.loocalCoordStart=self.seqChromEnd-self.mutChromEnd
			self.loocalCoordEnd=self.seqChromEnd-self.mutChromStart
		return self			
	
	def changeNtsBasedOnStrand(self):
		''' SNP on different strand, get revereComplement for wt, mut nt'''
		if self.sameStrand==False:
			self.refNts=reverseComp(self.refNts)
			self.altNts=reverseComp(self.altNts)
		return self

	def allPossiblePointMuts(self):
		''' Generate all possible point mutations for the sequence '''
		mutations = []
		allNts='ACGU'
		position = 0
		while (position < len(self.sequence)):
			nt = self.sequence[position]
			NtChoices = allNts.replace(nt,"")
			for NtChoice in NtChoices:
				mutation = nt+str(position+1)+NtChoice
				mutations.append(mutation)
			position+=1
		return mutations
	
	def randomMutsSampling(self,sampling=0.10):
		''' Generate all possible point mutations, sample 
		a certain percentage '''
		AllPointMutations = self.allPossiblePointMuts()
		SampleSize = int(sampling * len(AllPointMutations))
		random.shuffle(AllPointMutations)
		PointMutationsSampling = AllPointMutations[:SampleSize]
		return PointMutationsSampling
		
		
	def SequenceSanityCheck(self):
		''' Make sure there no non-standard characters in sequence '''
		seqIndex = 0
		while(seqIndex<len(self.sequence)):
			if (self.sequence[seqIndex] not in self.ntsAllowed):
				print self.sequence[seqIndex]+" not allowed as a nucleotide. exiting..."
				sys.exit()
						
			seqIndex+=1
		return self

	def buildUpAln(self,wt,mutpos,mut):
		''' adjust current alignment of wt, mut seqs accordingly '''
		#pos = pos + self.insOffset
		#print "INPUT POSITION : ",str(mutpos)
		#print "STUFF"
		#print self.pairwiseAln
		#print self.InsPositions
		if self.mutation!=None:
			if self.pairwiseAln==None:
				self.pairwiseAln = [list(self.sequence),list(xrange(len(self.sequence))),list(self.sequence),list(xrange(len(self.sequence)))]
				#print self.pairwiseAln
			if mut=="-":
				self.pairwiseAln[2] = self.pairwiseAln[2][:mutpos+self.insOffset]+["-"]*len(wt)+self.pairwiseAln[2][self.insOffset+mutpos+len(wt):]
				self.pairwiseAln[3] = self.pairwiseAln[3][:mutpos+self.insOffset]+[-1]*len(wt)+self.pairwiseAln[3][self.insOffset+mutpos+len(wt):]
				for val in xrange(self.insOffset+mutpos+len(wt),len(self.pairwiseAln[3])):
					if self.pairwiseAln[3][val]!=-1:
						self.pairwiseAln[3][val]-=len(wt)
			elif wt=="-":
				self.pairwiseAln[0] = self.pairwiseAln[0][:mutpos+self.insOffset]+list("-"*len(mut))+self.pairwiseAln[0][self.insOffset+mutpos:]
				self.pairwiseAln[1] = self.pairwiseAln[1][:mutpos+self.insOffset]+[-1]*len(mut)+self.pairwiseAln[1][self.insOffset+mutpos:]
				self.pairwiseAln[2] = self.pairwiseAln[2][:mutpos+self.insOffset]+list(mut)+self.pairwiseAln[2][mutpos+self.insOffset:]
				for count in xrange(mutpos,len(self.pairwiseAln[3])):
					if self.pairwiseAln[3][count]!=-1:
						self.pairwiseAln[3][count]+=len(mut)
				self.pairwiseAln[3] = self.pairwiseAln[3][:mutpos+self.insOffset]+[-1]*len(mut)+self.pairwiseAln[3][mutpos+self.insOffset:]
			elif (wt!="-" and  mut!="-" and (len(wt)!=1 or len(mut)!=1)):
				wtExt,mutExt = sizeNormalizeSeqs(wt,mut)
				self.pairwiseAln[0]=self.pairwiseAln[0][:mutpos+self.insOffset]+list(wtExt)+self.pairwiseAln[0][mutpos+self.insOffset+len(wt):]
				self.pairwiseAln[2]=self.pairwiseAln[2][:mutpos+self.insOffset]+list(mutExt)+self.pairwiseAln[2][mutpos+self.insOffset+len(wt):]

				self.pairwiseAln[1][self.insOffset+mutpos+len(wtExt):] = [val - wtExt.count("-") for val in self.pairwiseAln[1][self.insOffset+mutpos+len(wtExt):]]
				self.pairwiseAln[1] = self.pairwiseAln[1][:mutpos+self.insOffset]+range(mutpos+self.insOffset,mutpos+self.insOffset+len(wt))+[-1]*(len(wtExt)-len(wt))+self.pairwiseAln[1][self.insOffset+mutpos+len(wtExt):]
				self.pairwiseAln[3][self.insOffset+mutpos+len(mutExt):] = [val - mutExt.count("-") for val in self.pairwiseAln[3][self.insOffset+mutpos+len(mutExt):]]
				self.pairwiseAln[3] = self.pairwiseAln[3][:mutpos+self.insOffset]+range(mutpos+self.insOffset,mutpos+self.insOffset+len(mut))+[-1]*(len(mutExt)-len(mut))+self.pairwiseAln[3][self.insOffset+mutpos+len(mutExt):]
				self.insOffset+=wtExt.count("-")
				
			else:
				self.pairwiseAln[2][mutpos+self.insOffset:mutpos+self.insOffset+len(mut)]= mut
		return self
			

	def mutate_sequence(self,wt,pos,mut):
		''' Switch wildtype Nt in sequence with mutant Nt '''
		if self.mutantSeq==None:
			self.mutantSeq=self.sequence
		if mut=="-":
			self.mutantSeq = self.mutantSeq[:pos]+"-"*len(wt)+self.mutantSeq[pos+len(wt):]
		elif wt=="-":
			self.mutantSeq = self.mutantSeq[:pos]+mut+self.mutantSeq[pos:]
		else:
			self.mutantSeq = self.mutantSeq[:pos]+mut+self.mutantSeq[pos+1:]
		return self
		

	def ProcessMutation(self):
		''' extract position, make sure position in sequence matches
		the listed wildtype nucleotide in mutation name '''
		import re
		SNPs = []
		

		mut_positions = dict()
		bad_SNPs=[]
		for SNP in self.mutation.split(","):
			error=False
			positionset=re.findall(r"\d+",SNP)
			if len(positionset)!=1:
				error = True
			else:
				mutpos=int(positionset[0])-1
				[wt,mut]=string.split(SNP,positionset[0])
				if wt == mut:
					error = True
				positions = []
				for ntset in (wt,mut):				
					offset = 0
					for nt in list(ntset):
						#print offset
						#print nt,mutpos+offset,positions
						if nt not in self.mutCharsAllowed:
							error=True
						elif mutpos+offset not in positions:
							#print mutpos+offset
							positions.append(mutpos+offset)
						offset+=1
				
				if (self.sequence[mutpos:mutpos+len(wt)]!=wt and wt!="-"):
					error=True
				
				for pos in positions:
					if pos in mut_positions:
						error = True
					else:
						mut_positions[pos]=1
				
				if wt == mut:
					error = True
				elif wt == "-":
					self.Insertions.append((wt,mutpos,mut))
				elif mut == "-":
					self.Deletions.append((wt,mutpos,mut))
				else:
					self.Substitutions.append((wt,mutpos,mut))

			if error == False:
				pass
			else:
				bad_SNPs.append(SNP)
	
		#quit()
		if len(bad_SNPs)>0:
			print "ERROR WITH THE FOLLOWING SNPs : "
			for item in bad_SNPs:
				print item
			sys.exit(1)
			self.badSNP=True
		import operator
		for mutset in (self.Substitutions,self.Deletions,self.Insertions):
			mutset.sort(key=operator.itemgetter(1),reverse=True)
			#print mutset
			for mutation in mutset:
				#print "MUTATION"
				#print mutation
				self.buildUpAln(mutation[0],mutation[1],mutation[2])


		
		self.SNPPositions.sort()
		self.InsPositions.sort()
		self.DelPositions = merge_intervals(self.DelPositions)	
		#print "ALN : "
		#print numpy.array([list(seq) for seq in self.pairwiseAln])
		self.mask_aln()
		return self
		

	def rfam_cmScan(self,EvalCutoff=10,EvalIncCutoff=None,cmDB=os.getenv('HOME')+"/RFAMdb_11.0/Rfam.cm.1_1",printLines=False):
		''' search rfam to see if sequence matches any rfam families in a .cm file '''
		outputData=[]
		sequence=get_interval(self.seqChrom,self.seqChromStart,self.seqChromEnd,strand=self.strand)
		tempIn = tempfile.NamedTemporaryFile()
		tempOut = tempfile.NamedTemporaryFile()
		file = open(tempIn.name,"w")
		if self.seqName!=None:
			name = self.seqName
		else:
			name = "SEQ"
		file.write(">"+name+"\n"+sequence+"\n")
		file.close()
		os.system(" ".join(["cmscan","-o","/dev/null","--tblout",tempOut.name,cmDB,tempIn.name]))		

		output = open(tempOut.name,"r")
		lines = output.readlines()
		#target name         accession query name                          accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
		for line in lines:
			if printLines==True:
				print line.rstrip()
			if line[0]=="#":
				continue
			lineData=string.split(line.rstrip())
			lineData[17]="_".join(lineData[17:])
			lineData=lineData[:18]
			outputData.append(lineData)

		tempIn.close()
		tempOut.close()

		return outputData
			  

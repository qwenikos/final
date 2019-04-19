##python ~/biothesis/github/TSS/6_validation/FINAL/9_RNASeq_CAGE_per_gene/1_Bam_short_index_clear.py --inFile ENCFF125YZE_genome_alignment.bam
import os.path
import argparse
import os
from math import *


########################################
def printHeader(note):
    print "***************************"
    print "*                         *"
    print "*  "+str(note)
    print "*                         *"
    print "***************************"
    return


#######################################
printHeader("clean Bam file from contigs ")


parser = argparse.ArgumentParser()
parser.add_argument("--inFile", help="input Bam Filename")
parser.add_argument("--outDir", help="output Sam Dir")

args = parser.parse_args()
if args.inFile:  #i if in File exitst
    inFile=args.inFile
    inFilename=os.path.basename(inFile)
    inputDir=os.path.dirname(inFile)+"/"
    if (inputDir=="/"):
        inputDir=""
else:
	print "please give input fileName"
	exit()

if args.outDir:  #i if in File exitst
    outDir=args.outDir
    outputDir=os.path.dirname(outDir)+"/"
    if (outputDir=="/"):
        outputDir=""
else:
    outputDir=inputDir

sortedInputFile=outputDir+"sorted_"+inFilename
systemCall="samtools sort  -@ 8 -m 24G -T .  -o "+sortedInputFile +" "+inFile
print systemCall
os.system(systemCall)

systemCall="samtools index "+sortedInputFile
print systemCall
os.system(systemCall)

outputSamFileName=outputDir+inFilename[:-4]+".sam"

systemCall="samtools view -H "+sortedInputFile+" > "+outputSamFileName
print systemCall
os.system(systemCall)

for region in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]:
    systemCall="samtools view "+sortedInputFile+" '"+region+"' >> "+outputSamFileName
    print systemCall
    os.system(systemCall)
    


##Delete intermediate files
systemCall="rm -f  "+sortedInputFile
print systemCall
os.system(systemCall)
systemCall="rm -f  "+sortedInputFile+".bai"
print systemCall
os.system(systemCall)


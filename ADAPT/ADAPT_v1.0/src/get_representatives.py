import os
import sys


if len(sys.argv)>2:
    bedFileName = sys.argv[1]
    representativeFileName= sys.argv[2]
    bedNewNameFileName= sys.argv[3]
else:
    print "Give all FileNames"
    exit()
 
cnt=0

bedFile=open(bedFileName,"r")
bedNewNameFile=open(bedNewNameFileName,"w")
representativeFile=open(representativeFileName,"w")
for aLine in bedFile:
    cnt+=1
    cols=aLine.rstrip().split("\t")
    start=int(cols[1])
    end=int(cols[2])
    center_start=start+((end-start)-(end-start)%2)/2
    center_end=center_start+1
    
    cols[3]="fantom_tc_"+str(cnt)
    newLine="\t".join(cols)
    bedNewNameFile.write(newLine+"\n")
    cols[1]=str(center_start)
    cols[2]=str(center_end)
    newLineRepr="\t".join(cols)
    representativeFile.write(newLineRepr+"\n")
bedNewNameFile.close()
bedFile.close()
representativeFile.close()
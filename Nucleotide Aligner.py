#!/usr/bin/env python
import numpy as np
import sys

filepath=sys.argv[1]

file=filepath
seq1Index=0
seq2Index=1
inGap=9
termGap=5
misMatch=6
match=0

f = open(file).read().splitlines()

#look through fasta file use ">" as identified of the sequence name, that index + 1 will be the sequence.
listindex = []
i=0

for i in range(len(f)):
    if f[i][:1] == ">":
        listindex.append(i)

#in case theres multiple sequences lets allow for a sequence index, save each sequence as a variable
seq1Index= seq1Index
seq2Index= seq2Index

seq1=f[listindex[seq1Index]+1]
seq2=f[listindex[seq2Index]+1]

#lets put gap values here so I don't mistype it somewhere.
termGap=termGap
inGap=inGap
misMatch=misMatch
match=match

#initialize the alignment scoring array, gap penalty at beginning of each seq is 5. Row 0 and column 0 assume all terminal gaps. rows=seq1, cols=seq2. array[col,row]
score = np.zeros((len(seq1)+1,len(seq2)+1))
score[0,0:len(seq2)+1]=np.arange(0,(len(seq2)+1)*termGap,termGap)
score[0:len(seq1)+1,0]=np.arange(0,(len(seq1)+1)*termGap,termGap)

#calculate the inside of the array. We want max matching values: 4 possible: Match, mismatch, gap in seq 1, gap in seq 2.

for j in range(1,score.shape[0]): #seq1 rows
    for k in range(1,score.shape[1]): #seq2 columns

        if k == score.shape[1]-1 and j == score.shape[0]-1:
            if seq1[j-1] == seq2[k-1]:
                a= score[j-1,k-1] + match
            else: a= score[j-1,k-1] + misMatch
            score[j,k]= min (score[j-1,k] + termGap, score[j,k-1] + termGap, a)
        #if match add match score
        elif seq1[j-1] == seq2[k-1]:
            score[j,k]=score[j-1,k-1] + match
        #elifs are for terminal gaps, if we're at max j or k
        elif j == score.shape[0]-1:
            score[j,k]= min (score[j-1,k] + inGap, score[j,k-1] + termGap, score[j-1,k-1] + misMatch)
        elif k == score.shape[1]-1:
            score[j,k]= min (score[j-1,k] + termGap, score[j,k-1] + inGap, score[j-1,k-1] + misMatch)
        else:
            score[j,k]= min (score[j-1,k] + inGap, score[j,k-1] + inGap, score[j-1,k-1] + misMatch)

#start traceback
seq1align=[]
seq2align=[]

#determine direction of final alignment point
j=score.shape[0]-1
k=score.shape[1]-1
tbp=len(seq1)+len(seq2)-1 #traceBackPosition -1 as index starts at 0

while j>0:
    #I can't use just minimum value, I need to make sure that I get the lowest value cell that's along the same "path" as the cell I'm in
    x = None
    z = None
    c = None

    if score[j-1,k-1] == score[j,k] - misMatch or score[j-1,k-1] == score[j,k] - match:
        z = score[j-1,k-1]

    if score[j,k-1] == score[j,k] - inGap or score[j,k-1] == score[j,k] - termGap:
        x = score[j,k-1]

    if score[j-1,k] == score[j,k] - inGap or score[j-1,k] == score[j,k] - termGap:
        c = score[j-1,k]
        
    l = [x,z,c]
    minscore = min(y for y in l if y is not None)

    if minscore == score[j,k]-match:
        seq1align.insert(0,seq1[j-1])
        seq2align.insert(0,seq2[k-1])
        tbp=tbp-1
        j=j-1
        k=k-1
    elif minscore == score[j,k]-misMatch:
        seq1align.insert(0,seq1[j-1])
        seq2align.insert(0,seq2[k-1])
        tbp=tbp-1
        j=j-1
        k=k-1
    elif minscore == score[j-1,k]:
        seq1align.insert(0,seq1[j-1])
        seq2align.insert(0,"-")
        tbp=tbp-1
        j=j-1
        k=k
    elif minscore == score[j,k-1]:
        seq1align.insert(0,"-")
        seq2align.insert(0,seq2[k-1])
        tbp=tbp-1
        j=j
        k=k-1
    else:
        if j != 0 or k != 0: print("Alignment error")
        break

print("Alignment score for",f[listindex[seq1Index]], "and", f[listindex[seq2Index]], "is", score[score.shape[0]-1,score.shape[1]-1])
print(f[listindex[seq1Index]],"".join(seq1align))
print(f[listindex[seq2Index]],"".join(seq2align))
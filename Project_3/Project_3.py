""" Name: Ali Eastman Oku
    Algorithmic Bioinformatics
    ZID: Z-1893417
    Project 3
"""

from __future__ import division
import sys


T = [0 for i in range(500000)]
tTotal = 0
cTotal = 0 #total computed aligned positions
correct = 0 #correctly aligned positions 
criteria = 5 #default critera value 
sensitivity = 0
specificity = 0


def usage_error():
    
    error = "Usage: python -m project3 [<fasta file> <MafFile>, optional: criteria <num>]" 
    
    if len(sys.argv) < 3:
        print("Pass two files")
        print(error)


def trueAlignment(seq1, seq2):
    global tTotal
    global T
    
    """
    Calculates the true alignment between two alignment strings 
    
    Parameters
    ----------
    seq1 : str 
        alignment sequence of the first species 
        
    seq2: str
        alignment sequence of the second species 
 
    """
    
    len1 = len(seq1)
    len2 = len(seq2)
    
    length = len1 #default length 
    
    if len2 < len1:
        length = len2
    
    tTotal = 0 #stores the total true aligned positions
    pos1 = 0 #current position in the first species 
    pos2 = 0 #current position in the second species 
    
    for i in range(length):
        if seq1[i].upper() in ["A", "C", "T", "G"]:   #checks that base is a valid base (not gaps)
            
            if seq2[i].upper() in ["A", "C", "T", "G"]:   #checks that base is a valid base (not gaps)
                
                T[pos1] = pos2  #pos2 of second species aligned to pos1 of first species 
                
            elif seq2[i] == "-": 
                T[pos1] = pos2 - 0.5 
            pos1 = pos1 + 1
            
            
        if seq2[i].upper() in ["A", "C", "T", "G"]:
            pos2 = pos2 + 1
            
    
    tTotal = pos1
    
    
def readTrueAlignmentFile (filename):
    
    """
    Reads in the true alignment file and calls a method that calculates the true 
    alignment between two  sequences 
    
    Parameters
    ----------
    filename : path to the file (fasta file) containing the sequences  
    
    Note
    ---------
    This functions calls the trueAlignment function 
 
    """
    
    print("Read and Process the FASTA file(True Alignment File)")
    
    
    #stores the the alignment information for the two species 
    #first position: species 1
    #second position: species 2 
    align_block = []
    
    count = 0
    
    okay = False
    with open(filename) as f:
        while True:
            okay = False
            try:
                line = next(f)

                if line[0] == ">":
                    line_1 = next(f)
                    align_block.append(line_1)
                    count = count + 1
                    
                okay = True

            except StopIteration:
                break  # End of file.          
    trueAlignment(align_block[0], align_block[1])
    

    
def computedAlignment(pos_1, pos_2, seq1, seq2):
    
    """
    Process the computed alignment from a given MAF file 
    
    Parameters
    ----------
    seq1 : str 
        alignment sequence of the first species 
        
    seq2: str
        alignment sequence of the second species 
 
    """
    
    global cTotal
    global correct 
    
    pos1 = pos_1
    pos2 = pos_2
    
    len1 = len(seq1)
    len2 = len(seq2)
    
    length = len1 #default length 
    
    if len2 < len1:
        length = len2
    
    for i in range(length):
        
        if pos1 >= tTotal:
            break
        
        if seq1[i].upper() in ["A", "C", "T", "G"]:
            y = T[pos1] #retrieve true aligned position informtion from array T
            
            if seq2[i].upper() in ["A", "C", "T", "G"]:
                y_1 = pos2 #pos2 of second species is aligned to pos1 of first species 
                
            elif seq2[i] == "-":
                y_1 = pos2 - 0.5 #pos2 is next base 
                
            diff = abs(y - y_1)
            
            if diff <= criteria:
                correct = correct + 1
                
            pos1 = pos1 + 1
            cTotal = cTotal + 1
            
        if seq2[i].upper() in ["A", "C", "T", "G"]:
            pos2 = pos2 + 1
            
            
def readMafFile (filename):
    
    """
    Reads in the computed alignment file (MAF file)  and calls a method that processes this alignment 
    between two sequences from the file 
    
    Parameters
    ----------
    filename : path to the file (maf file) containing the sequences  
    
    Note
    ---------
    This functions calls the computedAlignment function 
 
    """
    
    print("Read and Process MAF file (computed) file")
    
    okay = False
    with open(filename) as f:
        while True:
            okay = False
            try:
                line_1 = next(f)

                if line_1[0] == "s":
                    pos1 = int(line_1.split()[2])
                    human = line_1.split()[-1].upper()
                    

                    line_2 = next(f)
                    pos2 = int(line_2.split()[2])
                    chimp = line_2.split()[-1].upper()

                    computedAlignment(pos1, pos2, human, chimp)
                    
                okay = True

            except StopIteration:
                break  # End of file.
                
                
def computeMetrics(cor, tTot, cTot):
    
    global sensitivity
    global specificity 
    
    if tTot > 0:
        sensitivity = cor/tTot
        
    if cTot > 0:
        specificity = cor/cTot
    
    
def main(fastaFile, mafFile):
    
    global criteria 
    
    if len(sys.argv) == 4:
        criteria = int(sys.argv[3])
              
    readTrueAlignmentFile(fastaFile)
    readMafFile(mafFile)
    computeMetrics(correct, tTotal, cTotal)
    
    print("Criteria: ", criteria)
    print("True Match Count: ", tTotal)
    print("Computed Match Count: ", cTotal)
    print("Correct: ", correct)
    print("Sensitivty: ", sensitivity)
    print("Specificity: ", specificity)
    

usage_error()
main(sys.argv[1], sys.argv[2])

    
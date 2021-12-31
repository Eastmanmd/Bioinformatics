""" Name: Ali Eastman Oku
    Algorithmic Bioinformatics
    ZID: Z-1893417
    Project 2
"""

from __future__ import division
import sys
import numpy as np
from biotite.sequence.phylo import upgma


rows = 119
cols = 119
matrix = [[0 for i in range(cols)] for j in range(rows)]  #stores the distance matrix between the gemones 

genome_names = [] #stores the name and position of the genomes for traceback
high_sub_rate = 0

def usage_error():
    
    error = "Usage: file, <genome/gene>, <gene name>, <start> <end>"
    
    if len(sys.argv) < 2:
        print(error)
        sys.exit(0)
        
        
        
def compare(string1, string2):
    
    """
    Compares the individual bases between two sequences of two different species
    Note: This method ignores gaps
    
    Parameters
    ----------
    string1 : str
        The dna sequence of first species 
    string2 : str
        The dna sequence of first species 
    """
    
    match_no = 0
    mismatch_no = 0

    for i, base in enumerate(string1):
        
        if (string1[i] not in ["A", "C", "T", "G"]):
            continue

        if (string2[i] not in ["A", "C", "T", "G"]):
            continue
                     
        if (string1[i] == string2[i]):
            match_no = match_no + 1
            
        elif (string1[i] != string2[i]):
            mismatch_no = mismatch_no + 1
                
    return match_no, mismatch_no


def compare_alt(string1, string2):
    
    """
    Compares the individual bases between two sequences of two different species
    Note: This method ignores gaps
    
    Parameters
    ----------
    string1 : str
        The dna sequence of first species 
    string2 : str
        The dna sequence of first species 
    """
    string_len = len(string1)
    
    if string_len > len(string2):
        string_len = len(string2)
    
    match_no = 0
    mismatch_no = 0

    for i in range(string_len):
        
        if (string1[i] not in ["A", "C", "T", "G"]):
            continue

        if (string2[i] not in ["A", "C", "T", "G"]):
            continue
                     
        if (string1[i] == string2[i]):
            match_no = match_no + 1
            
        elif (string1[i] != string2[i]):
            mismatch_no = mismatch_no + 1
                
    return match_no, mismatch_no



def compute_substitution_rate(match, mismatch):
    
    global high_sub_rate
    

    """
    Computes the substitution rate between two species 
    
    Parameters
    ----------
    match : str
        The total count of base pair matches between species 

    mismatch : str
        The total count of base pair mismatches between species 
        
    returns: int
        The substitution rate
        
    """
    
    if match == 0:

        high_sub_rate = high_sub_rate + 1
        return 1
    
    substitution_rate = mismatch / (match + mismatch)
    
    return substitution_rate



def read_and_compare(filename1):
    
    """
    Reads each line in the file and extracts the lines corresponding to dna sequences 
    compares these sequences using the compare method. 
    
    Parameters
    ----------
    filename : str
        The maf file containing the sequences 
    """
    
    global matrix
    global genome_names 
    
    with open(filename1) as f1:
        lines = [line.rstrip() for line in f1]

    count_row = 0 #keeps track of the rows in the matrix
    for i, sequence in enumerate(lines):
        
        if sequence[0] == "s":
            
            seq1 = lines[i].split()[-1].upper()
            seqname = lines[i].split()[1].upper()
            
            genome_names.append(seqname)
        
            count_col = 0 #keeps tracks of the columns in the matrix 
            for j, sequence in enumerate(lines):
                
                if sequence[0] == "s":
            
                    seq2 = lines[j].split()[-1].upper()
            
                    match, mismatch = compare(seq1, seq2) 
                    substi_ratio = compute_substitution_rate(match, mismatch)
                    
                    matrix[count_row][count_col] = substi_ratio
                    
                    count_col = count_col + 1
                    
            count_row = count_row + 1
                    

            
def get_gene_sequence (name, start, end, sequence): 
    

    """
    Gets the gene sequence given the start and end positions of the gene 
    
    Parameters
    ----------
    name : str
        The name of the gene 

    start : int
        The start position of the gene 
        
    end: int
        The end position of the gene
        
    sequence: str
        The whole genoic sequence including gaps from the global alignment 
        
    """
    
    new_sequence = "" 
    #stores the updated sequence without the gaps 
    
    for base in sequence:
        
        if (base.isalpha() == True):  #gets rid of gaps
            
            new_sequence = new_sequence + base 
            
        else:
            continue 
            
    new_sequence = new_sequence[start-1:end] #since MAF files start at index 0 while 
    
    return new_sequence           


def display(distance_matrix):
    print("\n The Distance Matrix is:")
                                                
    for i in distance_matrix:
        for j in i:
            print("{:.13f}".format(j), end=" ")
        print()
        
        
def read_and_compare(filename1):
    
    """
    Reads each line in the file and extracts the lines corresponding to dna sequences 
    compares these sequences using the compare method. 
    
    Parameters
    ----------
    filename : str
        The maf file containing the sequences 
    """
    
    global matrix
    global genome_names 
    
    with open(filename1) as f1:
        lines = [line.rstrip() for line in f1]

    count_row = 0 #keeps track of the rows in the matrix
    for i, sequence in enumerate(lines):
        
        if sequence[0] == "s":
            
            seq1 = lines[i].split()[-1].upper()
            seqname = lines[i].split()[1].upper()
            
            genome_names.append(seqname)
        
            count_col = 0 #keeps tracks of the columns in the matrix 
            for j, sequence in enumerate(lines):
                
                if sequence[0] == "s":
            
                    seq2 = lines[j].split()[-1].upper()
            
                    match, mismatch = compare(seq1, seq2) 
                    substi_ratio = compute_substitution_rate(match, mismatch)
                    
                    matrix[count_row][count_col] = substi_ratio
                    
                    count_col = count_col + 1
                    
            count_row = count_row + 1
                                     
            
            
def main():
    
    if sys.argv[2] == "genome":
        read_and_compare(sys.argv[1])
        display(matrix)
        
    if sys.argv[2] == "gene":
        
        name = sys.argv[3]
        start = int(sys.argv[4])
        end = int(sys.argv[5])
        file = sys.argv[1]
        
        read_and_compare_gene(name, start, end, file)
        display(matrix)
        
        
usage_error()
main()
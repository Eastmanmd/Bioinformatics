""" Name: Ali Eastman Oku
    Algorithmic Bioinformatics
    ZID: Z-1893417
    Project 1
"""

from __future__ import division
import sys


rows = 4
cols = 4
matrix = [[0 for i in range(cols)] for j in range(rows)]

mismatch_no = 0
match_no = 0
transitions = 0
transversions = 0

def usage_error():
    
    error = "Usage: Pass the file"
    
    if len(sys.argv) < 2:
        print(error)
        sys.exit(0)

def compare(string1, string2):
    """
    Compares the individual bases between two sequences of two different species
    Note: This method ignores gaps and non-coding bases
    Parameters
    ----------
    string1 : str
        The dna sequence of first species 
    string2 : str
        The dna sequence of first species 
    """
        
    global match_no
    global mismatch_no
    global transitions
    global transversions 
    for i, base in enumerate(string1):
        
        if (string1[i] not in ["A", "C", "T", "G"]):
            continue

        if (string2[i] not in ["A", "C", "T", "G"]):
            continue
                     
        if string1[i] == "A":

            if string2[i] == "A":
                matrix[0][0] = matrix[0][0] + 1
                match_no = match_no + 1

            elif string2[i] == "C":
                matrix[0][1] = matrix[0][1] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1

            elif string2[i] == "T":
                matrix[0][2] = matrix[0][2] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1

            elif string2[i] == "G":
                matrix[0][3] = matrix[0][3] + 1
                mismatch_no = mismatch_no + 1
                transitions = transitions + 1


        elif (string1[i] == "C"):

            if (string2[i] == "A"):
                matrix[1][0] = matrix[1][0] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1

            elif (string2[i] == "C"):
                matrix[1][1] = matrix[1][1] + 1
                match_no = match_no + 1

            elif (string2[i] == "T"):
                matrix[1][2] = matrix[1][2] + 1
                mismatch_no = mismatch_no + 1
                transitions = transitions + 1

            elif (string2[i] == "G"):
                matrix[1][3] = matrix[1][3] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1


        elif (string1[i] == "T"):

            if (string2[i] == "A"):
                matrix[2][0] = matrix[2][0] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1

            elif (string2[i] == "C"):
                matrix[2][1] = matrix[2][1] + 1
                mismatch_no = mismatch_no + 1
                transitions = transitions + 1

            elif (string2[i] == "T"):
                matrix[2][2] = matrix[2][2] + 1
                match_no = match_no + 1

            elif (string2[i] == "G"):
                matrix[2][3] = matrix[2][3] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1


        elif (string1[i] == "G"):

            if (string2[i] == "A"):
                matrix[3][0] = matrix[3][0] + 1
                mismatch_no = mismatch_no + 1
                transitions = transitions + 1

            elif (string2[i] == "C"):
                matrix[3][1] = matrix[3][1] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1

            elif (string2[i] == "T"):
                matrix[3][2] = matrix[3][2] + 1
                mismatch_no = mismatch_no + 1
                transversions = transversions + 1

            elif (string2[i] == "G"):
                matrix[3][3] = matrix[3][3] + 1
                match_no = match_no + 1


def read_and_compare(filename):
    """
    Reads each line in the file and extracts the lines corresponding to dna sequences 
    compares these sequences using the compare method. 
    
    Parameters
    ----------
    filename : str
        The maf file containing the sequences 
    """
    okay = False
    with open(filename) as f:
        while True:
            okay = False
            try:
                line_1 = next(f)

                if line_1[0] == "s":
                    human = line_1.split()[-1].upper()

                    line_2 = next(f)
                    chimp = chimp = line_2.split()[-1].upper()

                    compare(chimp,human)
                okay = True

            except StopIteration:
                break  # End of file.
                
                
def compute_substitution_rate(match, mismatch):
    
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
    
    substitution_rate = mismatch / (match + mismatch)
    
    return substitution_rate


def transitions_to_transversion(transi_no, transver_no):

    """
    Computes transitions to transversion ratio between two species 
    
    Parameters
    ----------
    transi_no : integer
        The total number of transitions

    transi_no : integer
        The total number of transversions 
        
    returns: int
        The transition to transversions ratio 
        
    """
    ratio = transi_no/ transver_no
    
    return ratio


def display(mat):
    """
    Adds the dna bases to display final results
    
    Parameters
    ----------
    mat: 2D array 
        
    """
    print("DNA Comparison Matrix" + '\n')
    
    mat_copy = mat.copy()
    
    mat_copy[0].insert(0, "A")
    mat_copy[1].insert(0, "C")
    mat_copy[2].insert(0, "T")
    mat_copy[3].insert(0, "G")
    
    mat_copy.insert(0, [" ", "A", "C", "T", "G"])
    
    print('\n'.join([''.join(['{:<9}'.format(item) for item in row]) 
      for row in mat_copy]))

    
    
def main(filename):
    
    read_and_compare(filename)
    substitution_rate = compute_substitution_rate(match_no, mismatch_no)
    ti_to_tr = transitions_to_transversion(transitions, transversions)
    
    print("Number of transitions: ", transitions)
    print("Number of transversions: ", transversions) 
    print("Substitution Rate: ", substitution_rate)
    print("ti/tv: ", ti_to_tr)
    
    display(matrix)
    
    

usage_error()
main(sys.argv[1])

          


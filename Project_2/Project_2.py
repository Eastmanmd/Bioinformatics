""" Name: Ali Eastman Oku
    Algorithmic Bioinformatics
    ZID: Z-1893417
    Project 2
"""

from __future__ import division
import sys

mismatch_no = 0
match_no = 0
total_gap_count = 0 #total gaps found in the entire sequence 
len_gap_counts = [0 for i in range(1000)]  #stores the counts of gaps based on the length 

def usage_error():
    
    error = "Usage: Pass the file"
    
    if len(sys.argv) < 2:
        print(error)
        sys.exit(0)
        
        
        
def count_gaps(string1):
    
    """
    Parses through a sequence and counts the number of gaps in the given string sequence
    
    Parameters
    ----------
    string1 : str 
        The dna sequence of first species 
 
    """
    
    global total_gap_count
    global len_gap_counts
    continue_gap = False
    gap_len = 0
    
    for i, base in enumerate(string1):
        
        # If gap is encoutered 
        if base == "-":
            
            if (continue_gap) == False:
                total_gap_count = total_gap_count + 1
                continue_gap = True
                gap_len = gap_len + 1
                
            elif (continue_gap) == True:
                gap_len = gap_len + 1
                   
        else:
            if (continue_gap) == True:
                len_gap_counts[gap_len] = (len_gap_counts[gap_len]) + 1
                gap_len = 0        #resets the gap length back to zero
                continue_gap = False
                
    if (continue_gap) == True:
        len_gap_counts[gap_len] = (len_gap_counts[gap_len]) + 1
        gap_len = 0 
        continue_gap = False
                


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
    
    global match_no
    global mismatch_no
    
    for i, base in enumerate(string1):
        
        if (string1[i].isalpha() == False) or (string2[i].isalpha() == False):
            continue 
            
        else:
            
            if string1[i] == string2[i]:  #Basepair match, increase match count
                match_no = match_no + 1
                
            elif string1[i] != string2[i]:  #Basepair mismatch, increase the msimatch count
                mismatch_no = mismatch_no + 1 
                
                
                
def compute_gap_rate(gap_no, match_no, mismatch_no):
    
    """
    Computes the total gap rate between two species given the sequence information
    
    Parameters
    ----------
    gap_no : integer 
        The total number of gaps present in the entire file
        
    match_no : integer 
        The total number of matches present in the sequence file
        
    mismatch_no : integer 
        The total number of matches present in the sequence file
        
    ----------
    return: the gap rate between the two species 
 
    """
    
    gap_rate = gap_no / (match_no + mismatch_no + gap_no)
    
    return gap_rate


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
                    count_gaps(human)
                    

                    line_2 = next(f)
                    chimp = chimp = line_2.split()[-1].upper()
                    count_gaps(chimp)

                    compare(chimp,human)
                okay = True

            except StopIteration:
                break  # End of file.

                
def display(len_gap_counts):
    
    """
    Displays the counts for each individual gap length encountered in the file
    
    Parameters
    ----------
    len_gap_counts : list/matrix 
        Contains the counts of the individual gap lengths in the sequences 
    """
    
    print('{:<13}'.format("Gap Length"), '{:<7}'.format("Count"), '{:>13}'.format("Gap Frequency"))
    for i, items in enumerate(len_gap_counts):
        if items != 0:
            print('{:<13}'.format(i), '{:<7}'.format(items), '{:<13}'.format(round(items/sum(len_gap_counts), 4)))


def main(filename):
    
    read_and_compare(filename)
    gap_rate = compute_gap_rate(total_gap_count, match_no, mismatch_no)
    
    print("Gap Rate: ", gap_rate)
    print("Total Gap Count: ", total_gap_count)
    print("Number of matches: ", match_no)
    print("Number of mismatches: ", mismatch_no)
    
    display(len_gap_counts)

    
usage_error()
main(sys.argv[1])
                
                

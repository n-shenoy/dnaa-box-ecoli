#!/usr/bin/env python

""" Python script for locating the ori of the E coli. genome """

__author__ = "Navami Shenoy"


from collections import defaultdict

def reverse_complement(string):
    #finds the reverse complement of a DNA strand
    reversed_string = string[::-1]  # reverse the string
    rev_comp_string = ''
    
    base_complement = {'A':'T', 'T':'A', 
                        'G':'C', 'C':'G'}
    # find the complement of the reversed string
    for base in reversed_string:
        rev_comp_string = rev_comp_string + base_complement[base]

    return rev_comp_string 



def min_skew(string):
    #finds the positions in the sequence where the skew is minimum
    skew_val = 0
    skew_list = [0]
    min_pos = []
    scores = {'A':0, 'T':0, 'G':1, 'C':-1} 
    
    for base in string:
        skew_val += scores[base]
        skew_list.append(skew_val)
    
    min_val = min(skew_list)
    min_val_pos = skew_list.index(min_val)
    min_pos.append(min_val_pos)

    for i in range(0, len(skew_list)):
        if min_val == skew_list[i] and min_val_pos != i:
            min_pos.append(i) 

    return min_pos, min_val #return positions and minimum skew value



def Hamming_distance(seq1, seq2):
    # compute the Hamming distance between two strings 
    # assumes the strings have the same length
    num_mismatches = 0
    for i in range(0, len(seq1)):
        if not seq1[i] == seq2[i]:
            num_mismatches += 1

    return num_mismatches



def neighbor(pattern, num_mismatches, kmers): 
    #recursively generates the set of all k-mers whose Hamming distance
    #from a given pattern does not exceed a certain number of mismatches
    if num_mismatches == 0:
        kmers.add(pattern)
    else:
        bases = ['A', 'C', 'G', 'T']
        for i in range(0, len(pattern)):
            for j in range(0, len(bases)):
                new_kmer = pattern[:i] + bases[j] + pattern[i+1:]
                if num_mismatches <= 1:
                    kmers.add(new_kmer)
                else:
                    neighbor(new_kmer, num_mismatches-1, kmers)
                    



def most_frequent_kmer(text, k, max_mismatches): # k = length of the kmer
    #finds the most frequent kmers and their reverse complements
    #that occur in a text with upto a certain number of mismatches
    all_kmers = defaultdict(int)

    for i in range(0, len(text) - k + 1):
        kmers = set() #stores the "neighborhood" of each kmer
        neighbor(text[i:i+k], max_mismatches, kmers)

        for kmer in kmers:
            all_kmers[kmer] += 1

    for p in all_kmers.keys():
        for i in range(0, len(text) - k + 1):
            num_mismatches = Hamming_distance(text[i:i+k], reverse_complement(p))
            if num_mismatches <= max_mismatches:
                all_kmers[p] += 1

    most_freq_kmers = []
    for p in all_kmers.keys():
        if all_kmers[p] == max(all_kmers.values()):
            most_freq_kmers.append(p)
            most_freq_kmers.append(reverse_complement(p))
    
    return set(most_freq_kmers)



def score_mismatches(text, k, max_mismatches): # k = kmer length
    #scores each kmer by its number of mismatches 
    #highest score assigned to kmer with the least mismatches
    freq_with_mismatches = defaultdict(int)

    for i in range(0, len(text) - k + 1):
        kmers = set() 
        neighbor(text[i:i+k], max_mismatches, kmers)

        for kmer in kmers:
            freq_with_mismatches[kmer] = max_mismatches

    for kmer in freq_with_mismatches.keys():
        for i in range(0, len(text) - k + 1):
            num_mismatches = Hamming_distance(text[i:i+k], reverse_complement(kmer))
            if num_mismatches <= max_mismatches:
                freq_with_mismatches[kmer] = freq_with_mismatches[kmer] - num_mismatches

    scores = {}
    most_freq_kmers = most_frequent_kmer(text, k, max_mismatches)
    for kmer in most_freq_kmers:
        scores[kmer] = freq_with_mismatches[kmer]
        
    return scores
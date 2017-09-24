# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Subeen Kim

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return ('T')
    elif nucleotide == 'C':
        return ('G')
    elif nucleotide == 'G':
        return ('C')
    elif nucleotide == 'T':
        return ('A')
    else:
        return ('-')

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    compdna = ''
    revcomp = ''
    for base in dna:
        compdna = compdna + get_complement(base)
        dnalist = list(compdna)
        dnalist.reverse()
        revcomp = ''.join(dnalist)
    return revcomp

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    segnumb = len(dna)//3
    index = 0
    segments = []

    while index <= segnumb:
        codon = dna[(3*index):(3*(index+1))]
        if  codon in ['TAA', 'TAG', 'TGA']:
            break
        else:
            segments.append(codon)
            index = index+1

    snippetDNA = ''.join(segments)

    return snippetDNA

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    segnumb = len(dna)//3
    index = 0
    oneframe = []
    while index <= segnumb:

        fragment = dna[(3*index):]

        if dna[(3*index):(3*(index+1))] == 'ATG':
            oneframe.append(rest_of_ORF(fragment))
            index = index + len(rest_of_ORF(fragment))//3

        else:
            index = index + 1

    return oneframe

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    orfs = []

    index = 0

    while index < 3:
        dnastring = ''.join(find_all_ORFs_oneframe(dna[index:]))
        orfs.append(dnastring)
        index = index + 1

    return orfs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    allorfs = []

    for originalorfs in find_all_ORFs(dna):
        if originalorfs == '':
            continue
        else:
            allorfs.append(originalorfs)

    for reverseorfs in find_all_ORFs(get_reverse_complement(dna)):
        if reverseorfs == '':
            continue
        else:
            allorfs.append(reverseorfs)

    return allorfs

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("")
    ''
    """
    reference = 0
    longest = ''

    for strand in find_all_ORFs_both_strands(dna):

        if len(strand) > reference:
            reference = len(strand)
            longest = strand

        else:
            continue

    return longest

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        >>> longest_ORF_noncoding("ATGATGATG",100)
        9
    """

    index = 0
    maxlength = 0

    while index < num_trials:
        index += 1

        shuffledna = shuffle_string(dna)
        length = len(longest_ORF(shuffledna))

        #print(shuffledna,longest_ORF(shuffledna),find_all_ORFs_both_strands(shuffledna),find_all_ORFs(shuffledna),length)

        if length > maxlength:
            maxlength = length

        else:
            continue

    return maxlength

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region.

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    dna_AA ={'ATT':'I','ATC':'I','ATA':'I','CTT':'L','CTC':'L','CTA':'L','CTG':'L','TTA':'L','TTG':'L','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TTT':'F','TTC':'F','ATG':'M','TGT':'C','TGC':'C','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GGT':'G','GGC':'G','GGA':'G','GGG':'G','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S','TAT':'Y','TAC':'Y','TGG':'W','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','CAT':'H','CAC':'H','GAA':'E','GAG':'E','GAT':'D','GAC':'D','AAA':'K','AAG':'K','CGC':'R','CGT':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R'}

    proteinlen = len(dna)//3
    index = 0
    protein = []

    while index <= proteinlen:

        codon = dna[(3*index):(3*(index+1))]

        if codon in dna_AA:
            protein.append(dna_AA[codon])
            index = index+1

        else:
            index = index+1
            continue

    proteinstring = ''.join(protein)

    return proteinstring

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    threshold = longest_ORF_noncoding(dna, 100)
    aminoacids = []

    for sequence in find_all_ORFs_both_strands(dna):
        if len(sequence) >= threshold:
            aminoacids.append(coding_strand_to_AA(sequence))
        else:
            continue

    return aminoacids



if __name__ == "__main__":
    import doctest
    # doctest.testmod(verbose=True)
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
    #doctest.run_docstring_examples(,)

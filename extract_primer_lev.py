#!/usr/bin/env python
__author__ = "Billy Fournier"
__license__ = "GPL"
__email__ = "billyfournier2000@yahoo.com"
__credits__ = "William Walters"



# python extract_primers_lev.py <mappingfile> <merged.fastq> <output_name> <tol>
""" last field is edit distance tolerance
"""


from sys import argv
from string import upper
from itertools import product
from optparse import OptionParser

from skbio.sequence import DNA
from skbio.parse.sequences import parse_fastq
from skbio.format.sequences import format_fastq_record




""" Creates all combinations of a primer sequence with a known variable base/s.

EXAMPLE:
A primer "ATCGN" has 4 variations at the last position so 4 bases can be
represented. Those primers would be:
"ATCGA"     "ATCGT"     "ATCGC"     "ATCGG"
"""

def gen_primer_list(mapping_primer):
    var = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': 'AG', 'Y': 'CT',
             'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
             'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    if type(mapping_primer) is list:
        poss = [var[c] for c in mapping_primer[0]]
    if type(mapping_primer) is str:
        poss = [var[c] for c in mapping_primer]
    return list(product(*poss))




def get_primers(mapping_file):
    """ Generates and returns all forward and reverse primer combinations
    from the mapping file.
    mapping_data:  the sequenced data's mapping file
    Will raise error if either the LinkerPrimerSequence or ReversePrimer fields
        are not present
    """
    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': 'AG', 'Y': 'CT',
             'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
             'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    # processing mapping_file into header and mapping_data
    with open(mapping_file) as mappingfile:
        header = mappingfile.readline().split()
        mapping_data = []
        mappingfile.readline() ### Could cause issues. ###
        for line in mappingfile:
            mapping_data.append(list(line.split()))


    if "LinkerPrimerSequence" in header:
        primer_ix = header.index("LinkerPrimerSequence") # ix = index
    else:
        raise IndexError(
            ("Mapping file is missing LinkerPrimerSequence field."))
    if "ReversePrimer" in header:
        rev_primer_ix = header.index("ReversePrimer")
    else:
        raise IndexError(("Mapping file is missing ReversePrimer field."))

    raw_forward_primers = set([])
    raw_forward_rc_primers = set([])
    raw_reverse_primers = set([])
    raw_reverse_rc_primers = set([])

    for line in mapping_data:
        if len(line) >= rev_primer_ix: #to prevent error with empty lines
            raw_forward_primers.update([upper(primer).strip() for
                                    primer in line[primer_ix].split(',')])
            raw_reverse_primers.update([str(DNA(primer).rc()) for
                                       primer in line[rev_primer_ix].split(',')])


    if not raw_forward_primers:
        raise ValueError(("No forward primers detected in mapping file."))
    if not raw_reverse_primers:
        raise ValueError(("No reverse primers detected in mapping file."))

# Finding the forward primers, or rc of reverse primers indicates forward
# read. Finding the reverse primer, or rc of the forward primers, indicates
# the reverse read, so these sets are merged.
    # raw_forward_primers.update(raw_reverse_rc_primers)
    # raw_reverse_primers.update(raw_forward_rc_primers)

    forward_primers = set([])
    reverse_primers = set([])
    # forward_rc_primers = set([])
    # reverse_rc_primers = set([])

    for sequence in raw_forward_primers:
        forward_primers |= set(gen_primer_list(sequence))
    for sequence in raw_reverse_primers:
        reverse_primers |= set(gen_primer_list(sequence))

    forward_primers = list(''.join(nucleo) for nucleo in forward_primers)
    reverse_primers = list(''.join(nucleo) for nucleo in reverse_primers)

    return forward_primers, reverse_primers

def editDistance(s1,s2):
    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
        distances = newDistances
    return distances[-1]

# print(editDistance("kitten","sitting"))
# print(editDistance("rosettacode","raisethysword"))

def editSearchForward(s1,s2,tolerance):
    windowSize = len(s1)
    start = 0
    good = [tolerance,'',-1]
    index = -1
    for i in range(30):
        if i-windowSize <= 0:
            start = 0
        else:
            start = i - windowSize
        ed = editDistance(s1,s2[start:i])
        # print s2[start:i]
        if ed < good[0]:
            good[0] = ed
            good[1] = [s2[start:i]]
            good[2] = i
            # print good
    return good[2]

def editSearchReverse(s1,s2,tolerance):
    windowSize = len(s1)
    seqLen = len(s2)
    start = 0
    good = [tolerance,'',-1]
    index = -1
    for i in xrange(seqLen, seqLen-24, -1):
        if i + windowSize >= seqLen:
            end = seqLen
        else:
            end = i + windowSize
        ed = editDistance(s1,s2[i:end])
        # print s2[start:i]
        if ed < good[0]:
            good[0] = ed
            good[1] = [s2[i:end]]
            good[2] = i
            # print good
    return good[2]

def remove_primers(input_fastq, output_fastq,for_primers,rev_primers):
    count = 0
    with open(input_fastq) as read, open(output_fastq, "w") as out_seqs:
        for label,seq,qual in parse_fastq(read):
            for primerF,primerR in zip(for_primers,rev_primers):
                start_slice = editSearchForward(primerF,seq,ed_tol)
                end_slice = editSearchReverse(primerR,seq,ed_tol)
                # print type(start_slice), '\t',end_slice

            if (start_slice != -1) and (end_slice != -1):
                print seq[start_slice:end_slice]
                curr_seq = seq[start_slice:end_slice]
                curr_qual = qual[start_slice:end_slice]
                formatted_fastq_line = format_fastq_record(label, curr_seq, curr_qual)
                out_seqs.write("%s" % (formatted_fastq_line))


mapping_file = argv[1]
merged_file = argv[2]
output_file = argv[3]
ed_tol = argv[4]


forward_primers, reverse_primers = get_primers(mapping_file)
remove_primers(merged_file,output_file,forward_primers,reverse_primers)

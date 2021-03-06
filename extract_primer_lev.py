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


import time
import functools
import collections
def lru_cache(maxsize = 255, timeout = None):
    class _LRU_Cache_class(object):
        def __init__(self, input_func, max_size, timeout):
            self._input_func        = input_func
            self._max_size          = max_size
            self._timeout           = timeout

            # This will store the cache for this function, format - {caller1 : [OrderedDict1, last_refresh_time1], caller2 : [OrderedDict2, last_refresh_time2]}.
            #   In case of an instance method - the caller is the instance, in case called from a regular function - the caller is None.
            self._caches_dict        = {}

        def cache_clear(self, caller = None):
            # Remove the cache for the caller, only if exists:
            if caller in self._caches_dict:
                del self._caches_dict[caller]
                self._caches_dict[caller] = [collections.OrderedDict(), time.time()]

        def __get__(self, obj, objtype):
            """ Called for instance methods """
            return_func = functools.partial(self._cache_wrapper, obj)
            return_func.cache_clear = functools.partial(self.cache_clear, obj)
            # Return the wrapped function and wraps it to maintain the docstring and the name of the original function:
            return functools.wraps(self._input_func)(return_func)

        def __call__(self, *args, **kwargs):
            """ Called for regular functions """
            return self._cache_wrapper(None, *args, **kwargs)
        # Set the cache_clear function in the __call__ operator:
        __call__.cache_clear = cache_clear


        def _cache_wrapper(self, caller, *args, **kwargs):
            # Create a unique key including the types (in order to differentiate between 1 and '1'):
            kwargs_key = "".join(map(lambda x : str(x) + str(type(kwargs[x])) + str(kwargs[x]), sorted(kwargs)))
            key = "".join(map(lambda x : str(type(x)) + str(x) , args)) + kwargs_key

            # Check if caller exists, if not create one:
            if caller not in self._caches_dict:
                self._caches_dict[caller] = [collections.OrderedDict(), time.time()]
            else:
                # Validate in case the refresh time has passed:
                if self._timeout != None:
                    if time.time() - self._caches_dict[caller][1] > self._timeout:
                        self.cache_clear(caller)

            # Check if the key exists, if so - return it:
            cur_caller_cache_dict = self._caches_dict[caller][0]
            if key in cur_caller_cache_dict:
                return cur_caller_cache_dict[key]

            # Validate we didn't exceed the max_size:
            if len(cur_caller_cache_dict) >= self._max_size:
                # Delete the first item in the dict:
                cur_caller_cache_dict.popitem(False)

            # Call the function and store the data in the cache (call it with the caller in case it's an instance function - Ternary condition):
            cur_caller_cache_dict[key] = self._input_func(caller, *args, **kwargs) if caller != None else self._input_func(*args, **kwargs)
            return cur_caller_cache_dict[key]


    # Return the decorator wrapping the class (also wraps the instance to maintain the docstring and the name of the original function):
    return (lambda input_func : functools.wraps(input_func)(_LRU_Cache_class(input_func, maxsize, timeout)))

# @lru_cache(maxsize=4095)
# @lru_cache(maxsize=8190) # may need to alter cashe size
@lru_cache(maxsize=16380)
def ld(s, t):
	if not s: return len(t)
	if not t: return len(s)
	if s[0] == t[0]: return ld(s[1:], t[1:])
	l1 = ld(s, t[1:])
	l2 = ld(s[1:], t)
	l3 = ld(s[1:], t[1:])
	return 1 + min(l1, l2, l3)

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




import cProfile
def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func



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


"""
@param s1   primer string
@param s2   sequence string
@param tolerance    maximum edit distances
@return index     Start or End index to use in sequence string
"""
# Search from the front. V2
def editSearchForward(s1,s2,tolerance):
    windowSize = len(s1)
    ed_check = int(tolerance)
    index = 0
    for i in range(int(tolerance)):
        # ed = editDistance(s1,s2[0: windowSize - i])
        ed = ld(s1,s2[0: windowSize - i])
        if ed < ed_check:
            ed_check = ed
            index = windowSize - i
    return index

# Search from back
def editSearchReverse(s1,s2,tolerance):
    windowSize = len(s1)
    seqEnd = len(s2) - windowSize
    ed_check = tolerance
    index = 0
    for i in range(int(tolerance)):
        # ed = editDistance(s1,s2[seqEnd + i: -1])
        ed = ld(s1,s2[seqEnd + i: -1])
        if ed < ed_check:
            ed_check = ed
            index = i + seqEnd
    return index


@do_cprofile
def remove_primers(input_fastq, output_fastq,for_primers,rev_primers, ed_tol):
    count = 0
    with open(input_fastq) as read, open(output_fastq, "w") as out_seqs:
        for label,seq,qual in parse_fastq(read):
            for primerF,primerR in zip(for_primers,rev_primers):
                start_slice = editSearchForward(primerF,seq,ed_tol)
                end_slice = editSearchReverse(primerR,seq,ed_tol)
                # print type(start_slice), '\t',end_slice

            if (start_slice != -1) and (end_slice != -1):
                curr_seq = seq[start_slice:end_slice]
                curr_qual = qual[start_slice:end_slice]
                formatted_fastq_line = format_fastq_record(label, curr_seq, curr_qual)
                out_seqs.write("%s" % (formatted_fastq_line))


mapping_file = argv[1]
merged_file = argv[2]
output_file = argv[3]
ed_tol = argv[4]


forward_primers, reverse_primers = get_primers(mapping_file)

remove_primers(merged_file,output_file,forward_primers,reverse_primers, ed_tol)

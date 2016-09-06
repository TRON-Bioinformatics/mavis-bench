import pysam
import itertools
from copy import copy
import re
from structural_variant.constants import *
import structural_variant.annotate as ann
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Alphabet import Gapped
from Bio import SeqIO
import numpy as np
from Bio.Data.IUPACData import ambiguous_dna_values
import warnings
import datetime
import time

def score_cigar(cigar_tuples):
    """
    for an input set of cigar tuples computes scoring minimums and returns a measure describing
    the alignment
    """
    putative_consecutive_matches = [[0, False]]
    total_events = 0

    for cigar, freq in cigar_tuples:
        if cigar == CIGAR.M:
            raise AttributeError('input cigar cannot have ambiguous values \'M\' is not allowed')
        if cigar != CIGAR.EQ and cigar != CIGAR.S:
            total_events += 1
        if cigar == CIGAR.EQ:
           putative_consecutive_matches[-1][0] += freq
        elif putative_consecutive_matches[-1][1]:
            pass
        else:
            putative_consecutive_matches.append([0, False])
    exact = [0]
    fuzzy = [0]
    for l, event in putative_consecutive_matches:
        if event:
            fuzzy.append(l)
        else:
            exact.append(l)
    return max(exact), max(fuzzy), total_events # longest exact, longest fuzzy, number of events

def longest_homopolymer(sequence):
    if len(sequence) == 0:
        raise AttributeError('cannot compute longest_homopolymer on an empty sequence')
    hp = [0]
    last_char = sequence[0]

    for char in sequence:
        if char == last_char:
            hp[-1] += 1
        else:
            hp.append(0)
        last_char = char
    return max(hp)

def str_cigar(read):
    """
    >>> read = pysam.AlignedSegment()
    >>> read.query_sequence = 'ATAGAACCCTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTAGAAATTCCACCACATGGAGGTTATTTGAG'
    >>> read.cigar = [(7, 55), (4, 70)]
    >>> str_cigar(read)
    '-----------------------------------------------------------------ATAGAACCCTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTagtgttccaaatgtgttggttctgcctagtatactgatttagaaattccaccacatggaggttatttgag--------------------------------------------------'
    """
    index = 0
    s = ''
    clipped = ''
    clipped_left = False
    max_str_len = 120
    primary_tuples = read.cigar

    typ, freq = read.cigar[0]
    typ_end, freq_end = read.cigar[-1]

    if (typ == CIGAR.S and typ_end != CIGAR.S) \
        or (typ == CIGAR.S and freq > freq_end): # soft clipped on the left side
        freq = read.cigar[0][1]
        clipped = read.query_sequence[0:freq]
        primary_tuples = read.cigar[1:]
        clipped_left = True
        index += freq
    elif typ_end == CIGAR.S: # soft clipped at the end
        freq = len(read.query_sequence) - read.cigar[-1][1]
        clipped = read.query_sequence[freq:]
        primary_tuples = read.cigar[:-1]
    
    for cigar, freq in primary_tuples:
        end = freq + index
        temp = read.query_sequence[index:end]
        if cigar == CIGAR.S:
            temp = temp.lower()
        index += freq
        s += temp
    
    if clipped_left:
        return clipped.rjust(max_str_len, '-').lower() + s.ljust(max_str_len, '-')
    return s.rjust(max_str_len, '-') + clipped.ljust(max_str_len, '-').lower()

def compute_cigar(ref, seq, **kwargs):
    """
    CIGAR VALUES
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch

    ASSUMES if force_end or force_start are given that all bases outside are SOFT CLIPPED
    # GTGAGTAAATTCAACATCGTTTTT
    # AACTTAGAATTCAAC---------
    >>> compute_cigar('GTGAGTAAATTCAACATCGTTTTT', 'AACTTAGAATTCAAC---------')
    [(4, 7), (7, 8)]
    >>> compute_cigar('TTTGCCATTGTTTGCTTTTGGAATCCAGGACACATCTTTAGGGCATGGCATTCCAGCTAACATCTTGGATGCCCCAGTGAAGATTTCAGTGCTTTTTTCCAATTTGTTTTTAGTATGCCTTGGCAAAAAAGAAAAGAGATACTAATTGATCAGAAAAAATGTTTTAGACTGCTAGCCCTGCTGTCTTTGGGGAAGTTGTATGCAGTGAGTAAATTCAACATCGTTTTTGGCCTCCCTATCAGTCATTAATGTAAAGTGGGGAGGCAGCTATTGCAGGCCACTATGATTTTGATAGTCAAGTAAAAGCTATGTTTTTTTGTTGCTGTTTGTTTATATCCATTAAGGGGAAAAATGGCCAGGCATGGTGGCTCACACCTGTAATCCCAGCACTTTGGGGGAGGCCAAGGCAGGAGGATCACTTGAGACCAGGAGTTTGAGACCACCCTGGGCAACATAGTGAGACCCTGTCTTTTTAAAGATAAATAAAGATTCTTTAAGAAAATAGAAAAGGAAAGGGGGAAAAATAATCCTCCTCACAGAATATTTGCAGTTCTTCTGTATGGAGAGAGGTAACCAAAACACATAATTATCTTCAGAATTTAGGATCATTGCTCTTTTTAAATTACATTCTATCCACGGGGATTTCTCAACAACTCAAAGGAGTAGTGTTGATTTATGGAGCTCAGCCCCTAGCAGTGTGCTAAAGCCCTAAGAAACCTGGCTCTGTTTTCATTTCTTGAGTTTAGCTCAAAAGAAGCTGCTTTTGGAGACATCTTAGGATATAGAACCCTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGAAGAGGTTTATTATAGGCTAACTATATAAGGTTATTTCAGCTGTGATGATCCAGCAGTTTTGTTGAACTTAAAAGAGCCTACCTATTAAGGATGCTTTATCGTGATGTAAAGAAGATGGTGCCTTGGTTAGTG', '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AGGGCTTTTTTTGTTGTTGTTATTT---------------------------------------------------------------------------------------------------------------------------------------')
    [(7, 25)]
    """
    if len(ref) != len(seq):
        raise AssertionError('input strings must be of equal length', ref, seq, len(ref), len(seq))
    assert(len(ref) > 0)
    
    force_start = kwargs.pop('force_start', None)
    force_end  = kwargs.pop('force_end', None)
    
    seq = Seq(str(seq), DNA_ALPHABET)
    ref = Seq(str(ref), DNA_ALPHABET)
    # ---XXXX gives you 3 as the start pos
    start = len(seq) - len(seq.lstrip(seq.alphabet.gap_char)) # don't need to add one b/c 0-indexed
    end = len(seq.rstrip(seq.alphabet.gap_char)) #----XXXX-- 10 - 2 = 8; exclusive 0-indexed range

    if force_start is not None and (force_start < start or force_start >= end):
        raise AttributeError('cannot force a start position unless it if >= to the aligned start position')
    if force_end is not None and (force_end > end or force_end <= start): 
        raise AttributeError('cannot force an end position past the length of the putative alignment')
    
    cigar_tuples = []
    if force_start is not None and force_start != start:
        cigar_tuples.append((CIGAR.S, force_start - start))
        start = force_start

    tail_sc = 0
    if force_end is not None and force_end != end:
        tail_sc = end - force_end
        end = force_end

    for i in range(start, end):
        current_mode = CIGAR.X
        if ref[i] == GAP and seq[i] == GAP:
            continue
        elif ref[i] == GAP:
            current_mode = CIGAR.I
        elif seq[i] == GAP:
            current_mode = CIGAR.D
        elif seq.alphabet.match(ref[i], seq[i]):
            current_mode = CIGAR.EQ
        
        if len(cigar_tuples) > 0:
            mode, freq = cigar_tuples[-1]
            if mode == current_mode:
                cigar_tuples[-1] = (current_mode, freq + 1)
            else:
                cigar_tuples.append((current_mode, 1))
        else:
            cigar_tuples.append((current_mode, 1))
    if tail_sc != 0:
        cigar_tuples.append((CIGAR.S, tail_sc))

    if force_start is None:
        # turn any X on the front into soft-clipping
        if cigar_tuples[0][0] == CIGAR.X:
            cigar_tuples[0] = (CIGAR.S, cigar_tuples[0][1])
    if force_end is None:
        if cigar_tuples[-1][0] == CIGAR.X:
            cigar_tuples[-1] = (CIGAR.S, cigar_tuples[-1][1])
    temp = sum([x[1] for x in cigar_tuples])
    if temp != len(seq.strip(GAP)):
        raise AssertionError('cigar tuple must be same length as input string', temp, len(ref), ref, cigar_tuples)
    return cigar_tuples

def build_string_from_reverse_path(ref, seq, path):
    """
    >>> ref = "-mxabdce"
    >>> seq = "-abc"
    >>> build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3)])
    ('e', '-')
    >>> build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3), (5,2), (4,2)])
    ('dce', '-c-')
    >>> build_string_from_reverse_path(ref, seq, [(6, 3), (5,2), (4,2), (3, 1), (2, 0)])
    ('mxabdce', '--ab-c-')
    """
    ref = str(ref)
    seq = str(seq)
    ref_string = []
    seq_string = []
    transitions = []

    end = (len(ref) -1, len(seq) - 1)
    if path[0] != end:
        path.insert(0, end)

    if path[-1][0] == 0 or path[-1][1] == 0:
        if path[-1] != (0, 0):
            path.append((0, 0))
    
    for i in range(0, len(path)):
        if i == 0:
            continue
        r, c = path[i]
        r0, c0 = path[i-1]
        up = r0 - r
        left = c0 - c
        if left == 0: # moving on ref and not seq
            seq_string.append('-'*up)
            temp = ref[r+1:r0+1]
            if up != len(temp):
                raise AssertionError('must be equal length', temp, up, path)
            ref_string.append(temp)
        elif up == 0:
            ref_string.append('-'*left)
            temp = seq[c+1:c0+1] 
            seq_string.append(temp)
            if left != len(temp):
                raise AssertionError('must be equal length', temp, left, path)
        else:
            if up != 1 or left != 1:
                raise AssertionError('diagonal moves must be up ({4}) one and left ({5}) one: {0},{1} => {2},{3}'.format(
                    r0, c0, r, c, up, left), path, i)
            ref_string.append(ref[r0])
            seq_string.append(seq[c0])
    ref_string.reverse()
    seq_string.reverse()
    return ''.join(ref_string), ''.join(seq_string)

def sw_pairwise_alignment(input_ref, input_seq, **kwargs):
    """
    aligns a short sequence to a relatively short reference sequence
    basic smith-waterman alignment for a pair of sequences (pairwise2 in biopython appears broken)

    arr[i, j] = max of the following
    - arr[i-1,j-1] + MATCH/MISMATCH PENALTY 
    - arr[i-1,j] + GAP_PENALTY
    - arr[i,j-1] + GAP_PENALTY

    deviates from the normal after computing the matrix, for all the best alignments recomputes their 
    scores ignoring consecutive gaps on either end (not penalized) but still penalizing interior gaps

    """
    max_inner_events = kwargs.pop('max_inner_events', 3)
    
    time_start = time.clock()
    size = (len(input_ref) + 1) * (len(input_seq) + 1) * 64 / 1000000
    
    if size > 10 and size < 2000:
        warnings.warn('warning the pairwise alignment matrix is large: {0} MB'.format(size))
    elif size > 2000:
        raise UserWarning('warning the pairwise alignment matrix is too large: {0} GB'.format(size/1000))
    
    ref = '-' + (input_ref.seq if hasattr(input_ref, 'seq') else input_ref)
    seq = '-' + (input_seq.seq if hasattr(input_seq, 'seq') else input_seq)
    
    arr = np.zeros(shape=(len(ref),len(seq)), dtype=np.int)

    MISMATCH = -1
    MATCH = 2
    GAP_PENALTY = -4
    highest_score = 0
    pointers = {} # key by location tuple
    # compute the dynamic programming matrix for the alignment
    for i in range(1, len(ref)):
        for j in range(1, len(seq)):
            pointers[(i, j)] = []
            matched = DNA_ALPHABET.match(ref[i], seq[j])
            diag = arr[i-1][j-1] + MATCH if matched else MISMATCH
            right = arr[i][j-1] + GAP_PENALTY
            down = arr[i-1][j] + GAP_PENALTY

            m = max([diag, right, down])
            arr[i][j] = m
            
            if diag == m:
                pointers[(i, j)].append((i-1, j-1))
            if right == m:
                pointers[(i, j)].append((i, j-1))
            if down == m:
                pointers[(i, j)].append((i-1, j))
            if j == len(seq) - 1:
                highest_score = max(arr[i][j], highest_score)

    # pull out the putative paths from the matrix
    paths = []
    for i in range(0, len(ref)): # must consume the smaller sequence
        j = len(seq) - 1
        if arr[i][j] == highest_score:
            paths.append([(i, j)])
    
    complete_paths = []
    while len(paths) > 0:
        current_path = paths.pop(0)
        s, t = current_path[-1] # last element of the path
        if s == 0 or t == 0:
            complete_paths.append(current_path)
            continue
        # check if we should continue with the current path based on the input
        if max_inner_events is not None and len(current_path) > max_inner_events + 1:
            # calculate the cigar string for the current path
            temp = build_string_from_reverse_path(ref, seq, current_path)
            c = compute_cigar(*temp)
            exact, fuzzy, events = score_cigar(c)
            if events > max_inner_events:
                continue
        for next_step in pointers[(s, t)]:
            paths.append(current_path + [next_step])

        # check
        if len(pointers[(s, t)]) == 0:
            print('ERROR: no back pointers')
            print('ref', ref)
            print('seq', seq)
            print('left', arr[s-1][t])
            print('up', arr[s][t-1])
            print('diag', arr[s-1][t-1])
            print('current', arr[s][t])
            print('coord', s, t, ref[s], seq[t])
            assert(len(next_step) != 0)
    
    # re-score the paths to give lower penalties to start and end gaps
    path_scores = []
    for path in complete_paths:
        ref_path = ''
        seq_path = ''
        #path.reverse()
        ref_path, seq_path = build_string_from_reverse_path(ref, seq, path)

        alignment_score = 0
        for i in range(0, len(ref_path)):
            if ref_path[i] == GAP or seq_path[i] == GAP:
                alignment_score += GAP_PENALTY
            else:
                alignment_score += MATCH if DNA_ALPHABET.match(ref_path[i], seq_path[i]) else MISMATCH
        temp = max([len(ref_path) - len(ref_path.lstrip(GAP)),
                    len(seq_path) - len(seq_path.lstrip(GAP))])
        alignment_score -= (temp - 1)*GAP_PENALTY
        temp = max([len(ref_path) - len(ref_path.rstrip(GAP)),
                    len(seq_path) - len(seq_path.rstrip(GAP))])
        alignment_score -= (temp - 1)*GAP_PENALTY
        path_scores.append((alignment_score, ref_path, seq_path, path))
    
    if len(path_scores) == 0:
        return []
    # choose the best alignment by score and then rightmost, and then sorted alphanumerically
    best_score = max([p[0] for p in path_scores])
    
    result = []
    for score, ref, seq, path in path_scores:
        if score != best_score:
            continue
        start = len(seq) - len(seq.lstrip(GAP)) # number of gaps, and 0-index
        
        a = pysam.AlignedSegment()
        a.query_sequence = seq.strip(GAP)
        a.reference_start = start 
        a.cigar = compute_cigar(ref, seq)
        result.append(a)
    t = (time.clock() - time_start)
    if t > 5:
        print('sw_pairwise_alignment:', t)
        print('# of paths investigated', len(path_scores))
        print('matrix size', (len(input_ref)+ 1) * (len(input_seq) + 1)/1000, 'K')
    return result

class Evidence:
    @property
    def break1(self):
        return self.breakpoint_pair.break1
    
    @property
    def break2(self):
        return self.breakpoint_pair.break2

    def _window(self, breakpoint):
        temp = self.read_length*2 + self.average_insert_size
        start = breakpoint.start - temp - self.call_error - self.read_length - 1
        end = breakpoint.end + temp + self.call_error + self.read_length - 1
        
        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + self.call_error + self.read_length - 1
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - self.call_error - self.read_length - 1
        return (start, end)
    
    def __init__(self, breakpoint_pair, bamfile, **kwargs):
        """
        @param =average_insert_size \a optional (type: int; default: 450)
        @param =min_anchor_size \a optional (type: int; default: 5)
        @param =min_mapping_quality \a optional (type: int; default: 20)
        @param =read_length \a optional (type: int; default: 125)
        """
        self.read_length = kwargs.pop('read_length', 125)
        self.average_insert_size = kwargs.pop('average_insert_size', 450)
        self.stdev_insert_size = kwargs.pop('stdev_insert_size', 25)
        self.call_error = kwargs.pop('call_error', 10)
        self.min_splits_reads_resolution = kwargs.pop('min_splits_reads_resolution', 3)
        self.convert_chr_to_index = kwargs.pop('convert_chr_to_index', {})
        self.bamfile = bamfile
        self.min_anchor_exact = kwargs.pop('min_anchor_exact', 6)
        self.min_anchor_fuzzy = kwargs.pop('min_anchor_fuzzy', 10) # allow a single event to interrupt the sequence
        self.max_anchor_events = kwargs.pop('max_anchor_events', 3)
        self.min_anchor_size = min(self.min_anchor_exact, self.min_anchor_fuzzy) # used when gathering evidence
        self.min_anchor_variation = kwargs.pop('min_anchor_variation', 0.3)
        self.min_mapping_quality = kwargs.pop('min_mapping_quality', 20)
        self.filter_secondary_alignments = kwargs.pop('filter_secondary_alignments', True)
        self.convert_index_to_chr = {}
        self.update_chr_to_index(self.convert_chr_to_index)

        self.breakpoint_pair = breakpoint_pair
        self.split_reads = {
                self.breakpoint_pair.break1: set(), 
                self.breakpoint_pair.break2: set()
                }
        self.flanking_reads = {
                self.breakpoint_pair.break1: set(), 
                self.breakpoint_pair.break2: set()
                }
    
    def add_flanking_read(self, read):
        if read.is_unmapped or read.mate_is_unmapped:
            raise UserWarning('input read (and its mate) must be mapped')
        w1 = self._window(self.break1)
        w1 = (w1[0] - 1, w1[1] - 1) # correct for psyam using 0-based coordinates
        w2 = self._window(self.break2)
        w2 = (w2[0] - 1, w2[1] - 1) # correct for psyam using 0-based coordinates
        
        # check if the strands are compatible
        if self.breakpoint_pair.opposing_strands is not None:
            opp = read.is_reverse == read.mate_is_reverse # same b/c natural orientation is opposite for paired reads
            if opp != self.breakpoint_pair.opposing_strands:
                raise UserWarning('strands are not compatible with expected strands')
        # check if this read falls in the first breakpoint window
        added_flanking = False
        if read.reference_start >= w1[0] and read.reference_end <= w1[1] \
                and read.reference_id == self.convert_chr_to_index[self.break1.chr]:
            if read.next_reference_start >= w2[0] and read.next_reference_start <= w2[1] \
                    and read.next_reference_id == self.convert_chr_to_index[self.break2.chr]:
                self.flanking_reads[self.break1].add(read)
                added_flanking = True
        if read.reference_start >= w2[0] and read.reference_end <= w2[1] \
                and self.convert_chr_to_index[self.break2.chr] == read.reference_id:
            if read.next_reference_start >= w1[0] and read.next_reference_start <= w1[1] \
                    and self.convert_chr_to_index[self.break1.chr] == read.next_reference_id:
                self.flanking_reads[self.break2].add(read)
                added_flanking = True
        
        if not added_flanking:
            raise UserWarning('does not map to the expected regions. does not support the current breakpoint pair')

    def add_split_read(self, read, first_breakpoint=True):
        
        breakpoint = self.break1 if first_breakpoint else self.break2
        opposite_breakpoint = self.break2 if first_breakpoint else self.break1
        
        if read.cigar[0][0] != CIGAR.S and read.cigar[-1][0] != CIGAR.S:
            return
        # the first breakpoint of a BreakpointPair is always the lower breakpoint
        # if this is being added to the second breakpoint then we'll need to check if the 
        # read soft-clipping needs to be adjusted
        
        if self.breakpoint_pair.stranded:
            if ( read.is_reverse and breakpoint.strand == STRAND.POS ) \
                    or ( not read.is_reverse and breakpoint.strand == STRAND.NEG ):
                        raise UserWarning('split read not on the appropriate strand')
        primary = ''
        clipped = ''
        if breakpoint.orient == ORIENT.LEFT:
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            clipped = read.query_sequence[read.query_alignment_end:] # end is exclusive in pysam
        elif breakpoint.orient == ORIENT.RIGHT:
            clipped = read.query_sequence[:read.query_alignment_start]
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
        else:
            raise AttributeError('cannot assign split reads to a breakpoint where the orientation has not been '
                    'specified')
        if len(primary) < self.min_anchor_size or len(clipped) < self.min_anchor_size:
            raise UserWarning('split read does not meet the minimum anchor criteria', primary, clipped)
        
        # try mapping the soft-clipped portion to the other breakpoint
        w = self._window(opposite_breakpoint)
        opposite_breakpoint_ref = ann.HUMAN_REFERENCE_GENOME[opposite_breakpoint.chr].seq[w[0] - 1: w[1]]
        a = sw_pairwise_alignment(opposite_breakpoint_ref, clipped)
        
        if len(a) != 1: # ignore multi-maps and poor alignments of clipped regions
            a = None
        else:
            a = a[0]
        # now that we've aligned the soft-clipped portion we should check if we need to
        # adjust the original read. only applies to the second breakpoint reads
        read = copy(read)
        # recalculate the read cigar string to ensure M is replaced with = or X
        read.cigar = self.recompute_cigar(read)
        if a is not None and not first_breakpoint:
            """
            then we need to ensure that the read preferentially aligns to the first breakpoint
            this matters for when the sequences at both breakpoints have common subsequence

            for example when the orientation of the original read is to the RIGHT
                    XX
                    --==========> (original read)
                    ||
                    ||
            =======>++            (re-aligned soft-clipped portion)

            or when the orientation is to the LEFT
                        XX
            ============->        (original read)
                        ||
                        ||
                        ++======> (re-aligned soft-clipped portion) 
            
            to account for this we take the matching bases from the aligned portion of
            the input read and mark them as soft-clipped, the read that is created from
            the soft-clipped portion gets these bases instead as aligned
            """
            shift = 0
            # loop through and check if the next base should be shifted
            while shift + a.reference_end <= len(opposite_breakpoint_ref) \
                    and shift + read.query_alignment_start < len(read.query_sequence) \
                    and shift + a.reference_start >= 0 \
                    and shift + read.query_alignment_end > read.query_alignment_start:
                
                new_ref_pos = a.reference_end + shift # first base after the alignment
                new_query_pos = shift + read.query_alignment_start # position in the primary sequence
                
                if breakpoint.orient == ORIENT.LEFT:
                    new_ref_pos = a.reference_start - 1 + shift # base before the alignment
                    new_query_pos = read.query_alignment_end - 1 + shift # last base of read alignment
                
                if not DNA_ALPHABET.match(opposite_breakpoint_ref[new_ref_pos], read.query_sequence[new_query_pos]): 
                    break # if they don't match we don't shift
                
                if shift + a.reference_end >= len(opposite_breakpoint_ref) \
                        or shift + read.query_alignment_start >= len(read.query_sequence)  \
                        or shift + a.reference_start <= 0 \
                        or shift + read.query_alignment_end <= read.query_alignment_start:
                            break # stop from over incrementing or over decrementing
                shift += 1 if breakpoint.orient == ORIENT.RIGHT else -1
            # at the end of the loop we have the number of bases we can 'shift' by
            if shift != 0:
                if breakpoint.orient == ORIENT.RIGHT:
                    # shift the original read start (and cigar) and the new read end (cigar)
                    clipped += primary[:shift]
                    primary = primary[shift:]
                    read.reference_start += shift
                    a.query_sequence = clipped
                    cigar = a.cigar[:] # can't update cigar by indexing, have to replace the whole list
                    if cigar[-1][0] == CIGAR.EQ:
                        cigar[-1] = (CIGAR.EQ, a.cigar[-1][1] + shift)
                    else:
                        cigar.append((CIGAR.EQ, shift))
                    a.cigar = cigar
                    read_reference = ann.HUMAN_REFERENCE_GENOME[self.convert_index_to_chr[read.reference_id]] \
                            .seq[read.reference_start - len(clipped):read.reference_end]
                    #read.cigar = compute_cigar(read_reference, read.query_sequence, force_start = len(clipped))
                    cigar = [(CIGAR.S, len(clipped))]
                    count = 0
                    for c, freq in read.cigar:
                        left = len(clipped) - count
                        if left > 0:
                            if freq > left:
                                cigar.append((c, freq - left))
                            count += freq
                        else:
                            cigar.append((c, freq))
                    read.cigar = cigar 
                else: # LEFT ====>------ to -----=====>
                    clipped = primary[len(primary) + shift:] + clipped # append to the front
                    primary = primary[:len(primary) + shift]
                    a.reference_start += shift
                    a.query_sequence = clipped
                    cigar = a.cigar[:] # can't update cigar by indexing, have to replace the whole list
                    if a.cigar[0][0] == CIGAR.EQ:
                        cigar[0] = (CIGAR.EQ, a.cigar[0][1] + abs(shift))
                    else:
                        cigar.insert(0, (CIGAR.EQ, abs(shift)))
                    a.cigar = cigar 
                    read_reference = ann.HUMAN_REFERENCE_GENOME[self.convert_index_to_chr[read.reference_id]] \
                            .seq[read.reference_start:read.reference_start + len(primary) + len(clipped)]
                    #read.cigar = compute_cigar(read_reference, read.query_sequence, force_end = len(primary))
                    cigar = [(CIGAR.S, len(clipped))]
                    count = 0
                    for c, freq in read.cigar[::-1]:
                        left = len(clipped) - count
                        if left > 0:
                            if freq > left:
                                cigar.insert(0, (c, freq - left))
                            count += freq
                        else:
                            cigar.insert(0, (c, freq))
                    read.cigar = cigar
                read.query_name += SUFFIX_DELIM + 'shift+{0}'.format(shift)
        if len(primary) < self.min_anchor_size or len(clipped) < self.min_anchor_size:
            raise UserWarning('split read does not meet the minimum anchor criteria')
        # update the read values
        if a is not None:
            if breakpoint.orient == ORIENT.LEFT: # right was clipped, left is on re-align
                a.query_sequence = primary + a.query_sequence
                a.cigar = [(CIGAR.S, len(primary))] + a.cigar
            else:
                a.query_sequence = a.query_sequence + primary
                a.cigar = a.cigar + [(CIGAR.S, len(primary))]
            w = self._window(opposite_breakpoint)
            a.reference_start = w[0] - 1 + a.reference_start  
            a.reference_id = self.convert_chr_to_index[opposite_breakpoint.chr]
            a.query_name = read.query_name + SUFFIX_DELIM + 'clipped-realign'
            a.next_reference_start = read.next_reference_start
            a.next_reference_id = read.next_reference_id
            a.flag = read.flag
            a.mapping_quality = NA_MAPPING_QUALITY # can't compute a mapping quality when we are only computing one alignment, 255=NA
        # need to do this after shifting
        s, t = self._window(breakpoint)
        s -= 1 # correct for pysam using 0-based coordinates
        t -= 1 # correct for pysam using 0-based coordinates
        if read.reference_start > t or read.reference_end < s \
                or self.convert_index_to_chr[read.reference_id] != breakpoint.chr:
            raise UserWarning('read does not map within the breakpoint evidence window')
        if len(read.query_sequence) - (read.query_alignment_end + 2) < self.min_anchor_size \
                and (read.query_alignment_start + 1) < self.min_anchor_size:
            raise UserWarning('split read does not meet the minimum anchor criteria')
        self.split_reads[breakpoint].add(read)
        if len(clipped) >= self.min_anchor_size and a is not None:
            self.split_reads[opposite_breakpoint].add(a)
    
    def recompute_cigar(self, read):
        """
        for cigar tuples where M is used, recompute to replace with X/= for increased
        utility and specificity

        CIGAR VALUES
        M 0 alignment match (can be a sequence match or mismatch)
        I 1 insertion to the reference
        D 2 deletion from the reference
        N 3 skipped region from the reference
        S 4 soft clipping (clipped sequences present in SEQ)
        H 5 hard clipping (clipped sequences NOT present in SEQ)
        P 6 padding (silent deletion from padded reference)
        = 7 sequence match
        X 8 sequence mismatch

        """
        temp = []
        offset = 0
        ref = ann.HUMAN_REFERENCE_GENOME[self.convert_index_to_chr[read.reference_id]].seq
        
        ref_pos = read.reference_start
        seq_pos = 0
        
        for cigar_value, freq in read.cigar:
            #print('current', cigar_value, freq)
            #print('ref_pos, seq_pos', ref_pos, seq_pos)
            #print(temp)
            if cigar_value in [CIGAR.S, CIGAR.I]:
                temp.append((cigar_value, freq))
                seq_pos += freq
            elif cigar_value in [CIGAR.D, CIGAR.N]:
                temp.append((cigar_value, freq))
                ref_pos += freq
            elif cigar_value in [CIGAR.M, CIGAR.X, CIGAR.EQ]:
                for offset in range(0, freq):
                    if DNA_ALPHABET.match(ref[ref_pos], read.query_sequence[seq_pos]):
                        if len(temp) == 0 or temp[-1][0] != CIGAR.EQ:
                            temp.append((CIGAR.EQ, 1))
                        else:
                            temp[-1] = (CIGAR.EQ, temp[-1][1] + 1)
                    else:
                        if len(temp) == 0 or temp[-1][0] != CIGAR.X:
                            temp.append((CIGAR.X, 1))
                        else:
                            temp[-1] = (CIGAR.X, temp[-1][1] + 1)
                    ref_pos += 1
                    seq_pos += 1
            else:
                raise NotImplementedError('unexpected CIGAR value {0} is not supported currently'.format(cigar_value))
        return temp
    

    def resolve_breakpoints(self):
        """
        use split read evidence to resolve bp-level calls for breakpoint pairs (where possible)
        if a bp level call is not possible for one of the breakpoints then returns None
        if no breakpoints can be resolved returns an empty list
        """
        pos1 = {}
        pos2 = {}
        
        for breakpoint, d in [(self.break1, pos1), (self.break2, pos2)]:
            for read in self.split_reads[breakpoint]:
                pos = breakpoint_pos(read, breakpoint.orient)
                # check if the read meets the minimum quality requirements
                exact, fuzzy, events = score_cigar(read.cigar)
                if (exact < self.min_anchor_exact and fuzzy < self.min_anchor_fuzzy) \
                        or events > self.max_anchor_events:
                            continue
                # check the variation in the aligned portion of the read
                algn = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
                hp = longest_homopolymer(algn)
                if 1 - hp / len(algn) < self.min_anchor_variation and hp > 3:
                    print('read anchor not variable', algn)
                    continue
                
                if pos not in d:
                    d[pos] = []
                d[pos].append(read)
            putative_positions = list(d.keys())
            for pos in putative_positions:
                if len(d[pos]) < self.min_splits_reads_resolution:
                    del d[pos]

        linked_pairings = []

        # now pair up the breakpoints with their putative partners
        for first, second in itertools.product(pos1, pos2):
            # can link a pair of breakpoints if the read prefix is in the evidence for both
            # or if the mate id is in the other breakpoint
            support = 0
            prefixes = {}
            for read in pos1[first]:
                prefix = read.query_name
                if SUFFIX_DELIM in read.query_name:
                    prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
                prefixes[prefix] = None
            for read in pos2[second]:
                prefix = read.query_name
                if SUFFIX_DELIM in read.query_name:
                    prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
                if prefix in prefixes:
                    support += 1
            if support == 0:
                continue
            linked_pairings.append((first, second, support, pos1[first], pos2[second]))
        
        for pos in pos1:
            if pos not in [ temp[0] for temp in linked_pairings ]:
                linked_pairings.append((pos, None, 0, pos1[pos], []))
        for pos in pos2:
            if pos not in [ temp[1] for temp in linked_pairings ]:
                linked_pairings.append((None, pos, 0, [], pos2[pos]))
        
        return linked_pairings
        
        """
        resolved_breakpoints = []
        for pos in seqs:
            if len(seqs[pos]) < self.min_splits_reads_resolution:
                continue
            mleft = max([len(l[0]) for l in seqs[pos]])
            mright = max([len(l[1]) for l in seqs[pos]])
            
            for i, s in enumerate(seqs[pos]):
                l, r, read = s
                l = str(l)
                r = str(r)
                seqs[pos][i] = SeqRecord(Seq(l.rjust(mleft, '-'), DNA_ALPHABET)) \
                        , SeqRecord(Seq(r.ljust(mright, '-'), DNA_ALPHABET)), read
            
            align = MultipleSeqAlignment([s[0] for s in seqs[pos]], alphabet=DNA_ALPHABET)
            summary = AlignInfo.SummaryInfo(align)
            lcons = summary.dumb_consensus(consensus_alpha = DNA_ALPHABET, ambiguous='N')
            align = MultipleSeqAlignment([s[1] for s in seqs[pos]], alphabet=DNA_ALPHABET)
            summary = AlignInfo.SummaryInfo(align)
            rcons = summary.dumb_consensus(consensus_alpha = DNA_ALPHABET, ambiguous='N')
            b = Breakpoint(
                    breakpoint.chr, 
                    pos, 
                    pos, 
                    strand = breakpoint.strand, 
                    orient = breakpoint.orient,
                    left_seq = lcons,
                    right_seq = rcons
                    ) 
            resolved_breakpoints.append((b, [x[2] for x in seqs[pos]]))
        return resolved_breakpoints"""
        return linked_pairings
    
    def update_chr_to_index(self, d):
        temp = {}
        for k, v in d.items():
            if v in temp:
                raise AttributeError('indicies must be unique', v)
            temp[v] = k
        self.convert_chr_to_index.update(d)
        self.convert_index_to_chr.update(temp)

    def load_evidence(self):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering
        """
         
        bamfile = pysam.AlignmentFile(self.bamfile, 'rb')
        convert_chr_to_index = {}
        for name in bamfile.references:
            convert_chr_to_index[name] = bamfile.gettid(name)
        
        self.update_chr_to_index(convert_chr_to_index)
        # TODO transcriptome window gathering
        
        for read in bamfile.fetch(
                '{0}'.format(self.break1.chr), 
                self._window(self.break1)[0],
                self._window(self.break1)[1]):
            if read.is_unmapped \
                    or read.mapping_quality < self.min_mapping_quality \
                    or (read.is_secondary and self.filter_secondary_alignments):
                continue
            try:
                self.add_split_read(read)
            except UserWarning:
                pass
            try:
                self.add_flanking_read(read)
            except UserWarning:
                pass
        for read in bamfile.fetch(
                '{0}'.format(self.break2.chr), 
                self._window(self.break2)[0],
                self._window(self.break2)[1]):
            if read.is_unmapped \
                    or read.mapping_quality < self.min_mapping_quality \
                    or (read.is_secondary and self.filter_secondary_alignments):
                continue
            try:
                self.add_split_read(read, False)
            except UserWarning:
                pass
            try:
                self.add_flanking_read(read)
            except UserWarning:
                pass

def breakpoint_pos(read, orient=ORIENT.NS):
    typ, freq = read.cigar[0]
    end_typ, end_freq = read.cigar[-1]
    ORIENT.enforce(orient)

    if typ != CIGAR.S and end_typ != CIGAR.S:
        raise AttributeError('cannot compute breakpoint for a read without soft-clipping')
    
    if orient == ORIENT.NS:
        if (typ == CIGAR.S and end_typ == CIGAR.S and freq > end_freq) \
            or typ == CIGAR.S and end_typ != CIGAR.S:
            orient = ORIENT.RIGHT
            # soft clipped to the left
        else:
            # soft clipped to the right
            orient = ORIENT.LEFT

    if orient == ORIENT.RIGHT:
        return read.reference_start
    else:
        return read.reference_end - 1

if __name__ == '__main__':
    import doctest
    doctest.testmod()

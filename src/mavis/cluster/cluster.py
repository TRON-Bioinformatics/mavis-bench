import itertools
from collections import namedtuple
from copy import copy
from typing import Dict, List

from ..breakpoint import Breakpoint, BreakpointPair
from ..constants import ORIENT, STRAND
from ..interval import Interval
from ..util import logger


class BreakpointPairGroupKey(
    namedtuple(
        'BreakpointPairGroupKey',
        [
            'chr1',
            'chr2',
            'orient1',
            'orient2',
            'strand1',
            'strand2',
            'opposing_strands',
            'explicit_strand',
        ],
    )
):
    def __new__(
        cls,
        chr1,
        chr2,
        orient1,
        orient2,
        strand1,
        strand2,
        opposing_strands=None,
        explicit_strand=False,
    ):
        if STRAND.NS in [strand1, strand2] and explicit_strand:
            raise ValueError('cannot have unspecified strand when explicit_strand is set')
        if not explicit_strand and opposing_strands is None:
            raise ValueError('opposing_strands must be specified when explicit_strand is false')
        if explicit_strand:
            opp = strand1 != strand2
            if opposing_strands is None:
                opposing_strands = opp
            elif opposing_strands != opp:
                raise ValueError(
                    'strand1 v strand2 v opposing_strands conflict.',
                    strand1,
                    strand2,
                    opposing_strands,
                )
        STRAND.enforce(strand1)
        STRAND.enforce(strand2)
        ORIENT.enforce(orient1)
        ORIENT.enforce(orient2)
        self = super(BreakpointPairGroupKey, cls).__new__(
            cls, chr1, chr2, orient1, orient2, strand1, strand2, opposing_strands, explicit_strand
        )
        return self


def weighted_mean(values, weights=None):
    if weights is None:
        weights = [1 for v in values]
    return sum(x * w for x, w in zip(values, weights)) / sum(weights)


def merge_integer_intervals(*intervals, weight_adjustment: int = 0) -> Interval:
    """
    Merges a set of integer intervals into a single interval where the center is the
    weighted mean of the input intervals. The weight is inversely proportional to the
    length of each interval. The length of the final interval is the average of the lengths
    of the input intervals capped in size so that it never extends beyond the union of the
    input intervals

    Args:
        weight_adjustment: add to length to lower weighting differences between small intervals
    """
    float_offset = 0.99999999
    intervals = list(intervals)
    centers = []
    weights = []
    lengths = []

    if not intervals:
        raise AttributeError(
            'cannot compute the weighted mean interval of an empty set of intervals'
        )
    for i in range(0, len(intervals)):
        curr = intervals[i]
        intervals[i] = Interval(curr[0], curr[1] + float_offset)
        for _ in range(0, intervals[i].freq):
            centers.append(intervals[i].center)
            weights.append((weight_adjustment + 1) / (intervals[i].length() + weight_adjustment))
            lengths.append(intervals[i].length())

    center = round(weighted_mean(centers, weights=weights) * 2, 0) / 2
    size = weighted_mean(lengths)  # -1 b/c center counts as one
    start = max([round(center - size / 2, 0), min([i[0] for i in intervals])])
    end = min([round(center + size / 2, 0), max([i[1] for i in intervals])])
    offset = min([center - start, end - center])
    result = Interval(
        int(round(center - offset, 0)), int(round(center + max(0, offset - float_offset), 0))
    )
    return result


def pair_key(pair):
    return (
        pair.break1.chr,
        pair.break2.chr,
        pair.break1.start,
        pair.break2.start,
        pair.break1.end,
        pair.break2.end,
        pair.break1.orient,
        pair.break2.orient,
        pair.break1.strand if pair.stranded else STRAND.NS,
        pair.break2.strand if pair.stranded else STRAND.NS,
        pair.stranded,
        pair.opposing_strands,
    )


def all_pair_group_keys(pair, explicit_strand=False):
    """
    Generate all possible group keys for a given breakpoint pair.

    This function generates all possible group keys for a given breakpoint pair. A group key represents a combination
    of various attributes of the breakpoint pair, such as chromosome, orientation, and strand.

    Args:
        pair (BreakpointPair): The breakpoint pair for which group keys are generated.
        explicit_strand (bool, optional): Flag indicating whether to consider explicit strand information.


    Returns:
        list: A list of BreakpointPairGroupKey objects representing all possible group keys.

    Example:
        >>> pair = BreakpointPair(...)
        >>> keys = all_pair_group_keys(pair)
        >>> print(keys)
        [BreakpointPairGroupKey(...), BreakpointPairGroupKey(...), ...]

    """
    opt = [
        [pair.break1.chr],
        [pair.break2.chr],
        ORIENT.expand(pair.break1.orient),
        ORIENT.expand(pair.break2.orient),
        [STRAND.NS] if not explicit_strand else STRAND.expand(pair.break1.strand),
        [STRAND.NS] if not explicit_strand else STRAND.expand(pair.break2.strand),
        [pair.opposing_strands],
    ]
    result = []
    for c1, c2, o1, o2, s1, s2, opp in list(itertools.product(*opt)):
        # Sanity check for the strand information
        if explicit_strand and (s1 != s2) != opp:
            continue
        elif not explicit_strand and opp == (o1 != o2):
            continue
        result.append(
            BreakpointPairGroupKey(c1, c2, o1, o2, s1, s2, opp, explicit_strand=explicit_strand)
        )
    return result


def merge_by_union(
    input_pairs: List[BreakpointPair],
    group_key: BreakpointPairGroupKey,
    weight_adjustment: int = 10,
    cluster_radius: int = 200,
) -> Dict[BreakpointPairGroupKey, List[BreakpointPair]]:
    """
    Merge the union of all breakpoint pairs that are within the given distance (cluster_radius).

    Args:
        input_pairs (List[BreakpointPair]): A list of breakpoint pairs to be merged.
        group_key (BreakpointPairGroupKey): The group key used for merging.
        weight_adjustment (int, optional): Weight adjustment factor. Defaults to 10.
        cluster_radius (int, optional): The maximum distance for merging pairs. Defaults to 200.

    Returns:
        Dict[BreakpointPairGroupKey, List[BreakpointPair]]: A dictionary where the keys are the group keys
        and the values are lists of merged breakpoint pairs.
    """
    # sort the pairs by start and end
    pairs_by_start = sorted(input_pairs, key=lambda x: x.break1.start)
    pairs_by_end = sorted(input_pairs, key=lambda x: x.break2.start)
    # create a dictionary of edges between pairs
    edges = {pair_key(p): set() for p in input_pairs}
    pairs_by_key = {}

    # create a dictionary of pairs by key
    for i in range(0, len(input_pairs)):
        pairs_by_key.setdefault(pair_key(pairs_by_start[i]), []).append(pairs_by_start[i])

        # To summarize, this line of code is populating a dictionary called pairs_by_key with key-value pairs. 
        # The keys are generated by calling the pair_key() function with pairs_by_start[i] as an argument, 
        # and the values are lists that contain the corresponding pairs_by_start[i] objects.
        
        # ordering contains the lists pairs_by_start and pairs_by_end after each other
        for ordering in [pairs_by_start, pairs_by_end]:
            # try all combinations until start distance alone is too far
            curr = ordering[i] # The first value of curr is the first element of pairs_by_start
            ckey = pair_key(curr)
            edges.setdefault(ckey, set())
            # Each breakpoint is compared to all other breakpoints in the list to determine if they are within the cluster_radius
            # The ordering process is just a improvement to the algorithm to reduce the number of comparisons
            for j in range(i + 1, len(input_pairs)):
                other = ordering[j]
                okey = pair_key(other)
                distance = abs(Interval.dist(curr.break1, other.break1))
                if distance > cluster_radius:
                    break
                distance += abs(Interval.dist(curr.break2, other.break2))
                if distance <= cluster_radius:
                    edges[ckey].add(okey)
                    edges[okey].add(ckey)

    # edges contains a key for each pair and a set of keys for all pairs within cluster_radius
    merged = set()
    merge_nodes = []
    for node in edges:              # node = pair_key(BreakpointPair) : set of pair_key(BreakpointPair)
        if node in merged:
            continue
        adj = edges[node] | {node}
        merged.add(node)
        unmerged = adj - merged
        # follow edges to merge all connected nodes until all edges have been visited
        # extracts the current connected component
        while unmerged:
            for other in unmerged:
                adj.update(edges[other])
                merged.add(other)
            unmerged = adj - merged
        merge_nodes.append(adj)
    
    nodes = {}
    for node_keys in merge_nodes:
        pairs = []
        for pkey in node_keys:
            pairs.extend(pairs_by_key[pkey])
        itvl1 = merge_integer_intervals(
            *[p.break1 for p in pairs], weight_adjustment=weight_adjustment
        )
        itvl2 = merge_integer_intervals(
            *[p.break2 for p in pairs], weight_adjustment=weight_adjustment
        )
        if group_key.chr1 == group_key.chr2:
            itvl1.end = min(itvl2.end, itvl1.end)
            itvl2.start = max(itvl2.start, itvl1.start)
            itvl1.start = min(itvl1.start, itvl1.end)
            itvl2.end = max(itvl2.end, itvl2.start)
        b1 = Breakpoint(
            group_key.chr1,
            itvl1.start,
            itvl1.end,
            orient=group_key.orient1,
            strand=group_key.strand1,
        )
        b2 = Breakpoint(
            group_key.chr2,
            itvl2.start,
            itvl2.end,
            orient=group_key.orient2,
            strand=group_key.strand2,
        )
        # create the new bpp representing the merge of the input pairs
        new_bpp = BreakpointPair(
            b1, b2, opposing_strands=group_key.opposing_strands, stranded=group_key.explicit_strand
        )
        nodes.setdefault(new_bpp, []).extend(pairs)
    return nodes


def merge_breakpoint_pairs(
    input_pairs: List[BreakpointPair],
    cluster_radius: int = 200,
    cluster_initial_size_limit: int = 25,
    verbose: bool = False,
) -> Dict[BreakpointPair, List[BreakpointPair]]:
    """
    Merges breakpoint pairs based on their proximity and size.

    This function performs a two-step merging process to merge breakpoint pairs. In the first step, it merges all 'small' events
    (defined by the cluster_initial_size_limit parameter) as the union of all events that fall within the cluster_radius parameter.
    For example, if two events are within the cluster_radius of each other, they are merged into a single event.

    In the second step, for all remaining events, the function chooses the 'best' merge for any event within the cluster_radius of an existing node.
    Otherwise, the node is added unmerged. The events in the second phase are processed in order of smallest total breakpoint interval size to largest.

    Args:
        input_pairs (List[BreakpointPair]): The breakpoint pairs to be merged.
        cluster_radius (int, optional): The maximum distance allowed for a node to merge. Defaults to 200.
        cluster_initial_size_limit (int, optional): The maximum size of breakpoint intervals allowed in the first merging phase. Defaults to 25.
        verbose (bool, optional): Whether to print verbose output during the merging process. Defaults to False.

    Returns:
        Dict[BreakpointPair, List[BreakpointPair]]: A mapping of merged breakpoint pairs to the input pairs used in the merge.
    """

    # Calculate the distance between the centers of two breakpoint pairs
    def pair_center_distance(pair1, pair2):
        d = abs(pair1.break1.center - pair2.break1.center)
        d += abs(pair1.break2.center - pair2.break2.center)
        return d

    mapping = {}
    groups = {}  # split the groups by putative pairings
    pair_weight = {}
    explicit_strand = False
    phase2_groups = {}
    for pair in input_pairs:
        if pair.stranded:
            explicit_strand = True
            break

    doubled = 0
    for i, old_pair in enumerate(input_pairs):
        pair = copy(old_pair)
        pair.data['tag'] = i
        # the pair_key function provides a way to generate a unique key for a pair object based on its attributes, 
        # allowing for efficient comparison and grouping of pairs in the project.
        k = pair_key(pair)

        # setdefault is used to retrieve the values of keys. Is a key already present in the dictionary, setdefault returns its corresponding value.
        # However, if a key is not present in the dictionary, setdefault adds it with the defined default value and returns that value.
        pair_weight.setdefault(k, []).append(pair)

        # BreakpointPairGroupKeys are generated for each BPP. It is a named tuple that represents a combination of various attributes
        # of the breakpoint pair, such as chromosome, orientation, and strand. It is used to compare and group pairs.
        putative_group_keys = all_pair_group_keys(pair, explicit_strand=explicit_strand)
        doubled += len(putative_group_keys)
        if len(putative_group_keys) < 1:
            raise NotImplementedError('bad breakpoint input does not fit any groups', pair)
        for key in putative_group_keys:
            # The length of the breakpoint is defined as the length of the interval of the breakpoint. The interval describes uncertainty also refered as confidence interval (CIPOS/CIEND).
            # If the sum of the lengths of the two breakpoints is greater than the cluster_initial_size_limit. This may be due to insertions or microhomologies.
            # The pair is added to the phase2_groups dictionary.
            if len(pair.break1) + len(pair.break2) > cluster_initial_size_limit:
                phase2_groups.setdefault(key, []).append(pair)
            else:
                groups.setdefault(key, []).append(pair)
    # now try all pairwise combinations within groups
    for group_key in sorted(set(list(groups) + list(phase2_groups))):
        # This calculates the number of pairs in each group for the given group_key.
        count = len(groups.get(group_key, [])) + len(phase2_groups.get(group_key, []))
        if verbose:
            logger.info(f'{group_key} pairs: {count}')
        nodes = merge_by_union(
            groups.get(group_key, []), # Returns the value of the specified key. If the key does not exist, an empty tuple.
            group_key,
            weight_adjustment=cluster_initial_size_limit,
            cluster_radius=cluster_radius,
        )

        # phase 2. Sort all the breakpoint pairs left by size and merge the smaller ones in first
        # this is be/c we assume that a larger breakpoint interval indicates less certainty in the call
        phase2_pairs = sorted(
            phase2_groups.get(group_key, []),
            key=lambda p: (len(p.break1) + len(p.break2), pair_key(p)),
        )

        for pair in phase2_pairs:
            # Berechne die Distanz zwischen dem Zentrum des aktuellen Breakpoint-Paares und allen anderen Breakpoint-Paaren in den zuvor erstellten nodes
            distances = sorted(
                [(pair_center_distance(pair, node), node) for node in nodes], key=lambda x: x[0]
            )
            merged = False

            if len(distances) > 0:
                # Bestimme den Knoten mit dem kleinsten Abstand, der kleiner ist als der cluster_radius
                best = min(distances, key=lambda x: x[0])
                for dist, node in distances:
                    if dist > best[0] or dist > cluster_radius:
                        break
                    # Verwende alle Knoten zum Mergen incl. die zuvor gemergten Knoten
                    pairs = nodes[node] + [pair]

                    itvl1 = merge_integer_intervals(
                        *[p.break1 for p in pairs], weight_adjustment=cluster_initial_size_limit
                    )
                    itvl2 = merge_integer_intervals(
                        *[p.break2 for p in pairs], weight_adjustment=cluster_initial_size_limit
                    )
                    if group_key.chr1 == group_key.chr2:
                        itvl1.end = min(itvl2.end, itvl1.end)
                        itvl1.start = min(itvl1.start, itvl1.end)
                        itvl2.start = max(
                            itvl2.start,
                            itvl1.start + 2 if not any([p.opposing_strands for p in pairs]) else 1,
                        )  # for merging putative deletion events
                        itvl1.start = min(itvl1.start, itvl1.end)
                        itvl2.end = max(itvl2.end, itvl2.start)

                    b1 = Breakpoint(
                        group_key.chr1,
                        itvl1.start,
                        itvl1.end,
                        orient=group_key.orient1,
                        strand=group_key.strand1,
                    )
                    b2 = Breakpoint(
                        group_key.chr2,
                        itvl2.start,
                        itvl2.end,
                        orient=group_key.orient2,
                        strand=group_key.strand2,
                    )

                    new_bpp = BreakpointPair(
                        b1,
                        b2,
                        opposing_strands=group_key.opposing_strands,
                        stranded=explicit_strand,
                    )
                    del nodes[node]
                    nodes.setdefault(new_bpp, []).extend(pairs)
                    merged = True
            if not merged:
                b1 = Breakpoint(
                    group_key.chr1,
                    pair.break1.start,
                    pair.break1.end,
                    orient=group_key.orient1,
                    strand=group_key.strand1,
                )

                b2 = Breakpoint(
                    group_key.chr2,
                    pair.break2.start,
                    pair.break2.end,
                    orient=group_key.orient2,
                    strand=group_key.strand2,
                )

                new_bpp = BreakpointPair(
                    b1, b2, opposing_strands=group_key.opposing_strands, stranded=explicit_strand
                )
                nodes.setdefault(new_bpp, []).append(pair)
        if verbose:
            logger.info(f'merged {count} down to {len(nodes)}')
        for node, pairs in nodes.items():
            if node in mapping:
                raise KeyError('duplicate merge node', str(node), node, pair_key(node))
            mapping[node] = pairs
    # assertion to check that no nodes were left out of merging
    merge_sources = set()
    for merge_node, sources in mapping.items():
        merge_sources.update([p.data['tag'] for p in sources])
    if len(merge_sources) != len(input_pairs):
        raise AssertionError(
            'merged node inputs ({}) does not equal the number of pairs input ({})'.format(
                len(merge_sources), len(input_pairs)
            )
        )
    return mapping

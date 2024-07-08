import collections
import dendropy
import math
import numpy


def computeDistancesOfTrees(reference_tree, comparison_tree, n_populations):
    reference_tree.encode_bipartitions()
    comparison_tree.encode_bipartitions()

    ref_edge_lengths = collections.defaultdict(float)
    cmp_edge_lengths = collections.defaultdict(float)

    combined_splits = set()

    reference_tree_length = 0.0

    for ref_edge in reference_tree.postorder_edge_iter():
        if ref_edge.head_node != reference_tree.seed_node: # ignore root edge
            ref_edge_lengths[ref_edge.bipartition.split_bitmask] += ref_edge.length
            combined_splits.add(ref_edge.bipartition.split_bitmask)
            reference_tree_length += ref_edge.length

    for cmp_edge in comparison_tree.postorder_edge_iter():
        if cmp_edge.head_node != comparison_tree.seed_node: # ignore root edge
            cmp_edge_lengths[cmp_edge.bipartition.split_bitmask] += cmp_edge.length
            combined_splits.add(cmp_edge.bipartition.split_bitmask)

    n_combined_splits = len(combined_splits)

    rf_components = numpy.zeros(n_combined_splits)
    weighted_rf_components = numpy.zeros(n_combined_splits)
    branch_score_components = numpy.zeros(n_combined_splits)

    i = 0

    for split_bitmask in combined_splits:
        if split_bitmask not in cmp_edge_lengths: # only in reference tree
            rf_components[i] += 1.0
            weighted_rf_components[i] += ref_edge_lengths[split_bitmask]
            branch_score_components[i] += ref_edge_lengths[split_bitmask]
        elif split_bitmask not in ref_edge_lengths: # only in comparison tree
            rf_components[i] += 1.0
            weighted_rf_components[i] += cmp_edge_lengths[split_bitmask]
            branch_score_components[i] += cmp_edge_lengths[split_bitmask]
        else: # in both
            branch_score_components[i] += abs(ref_edge_lengths[split_bitmask] - cmp_edge_lengths[split_bitmask])

        i += 1

    rf_distance = numpy.sum(rf_components)
    normalized_rf_distance = rf_distance / float((n_populations - 3) * 2)
    weighted_rf_distance = numpy.sum(weighted_rf_components)
    euclidean_branch_score = numpy.sqrt(numpy.sum(branch_score_components * branch_score_components))
    normalized_branch_score = numpy.sum(branch_score_components) / reference_tree_length


    print(
        rf_distance,
        normalized_rf_distance,
        weighted_rf_distance,
        euclidean_branch_score,
        normalized_branch_score)

    return rf_distance, normalized_rf_distance, weighted_rf_distance, euclidean_branch_score, normalized_branch_score

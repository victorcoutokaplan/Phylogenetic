import collections
import dendropy
import math
import numpy


def makeDendroTree(birth_rate, death_rate, n_populations, extra_populations):
    """
    Given birth_rate, death_rate, n_populations, extra_populations arguments,
    simulates a tree made under those params and returns it as a dendropy tree

    Returns: dendropy tree
    """

    # Make tree

    reference_tree = dendropy.simulate.treesim.birth_death_tree(
        birth_rate = birth_rate,
        death_rate = death_rate,
        gsa_ntax = n_populations + extra_populations,
        num_extant_tips = n_populations)


    reference_tree.is_rooted = True # seems to be important for compatability with msprime
    reference_tree.encode_bipartitions()
    # ref_edge_lengths = collections.defaultdict(float)
    # combined_splits = set()
    # reference_tree_length = 0.0
    #
    # # adding edge lengths
    #
    # for ref_edge in reference_tree.postorder_edge_iter():
    #     if ref_edge.head_node != reference_tree.seed_node: # ignore root edge
    #         ref_edge_lengths[ref_edge.bipartition.split_bitmask] += ref_edge.length
    #         combined_splits.add(ref_edge.bipartition.split_bitmask)
    #         reference_tree_length += ref_edge.length

    return reference_tree

import numpy as np
from typing import List
from scipy.stats import mannwhitneyu
from sklearn.model_selection import train_test_split


def create_jitter(x, noise=1e-7):
    "create jitter in the case of 0s. this would satisfy the Matrix/0 LinAlg error"
    jitter = np.ones(len(x))
    halfway = round(len(x) * 0.5)
    jitter[:halfway] = noise
    jitter[halfway:] = -1 * noise
    np.random.shuffle(jitter)
    return jitter


def ranges(nums: list) -> List[tuple]:
    """Given list of integers, returns list with consecutive integers as ranges.
    If list has single number range `x`, will return `(x, x)`.
    For Input: [2, 3, 4, 7, 8, 9, 15]
    Output looks like: [(2, 4), (7, 9), (15, 15)]
    """
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))


def find_intersecting_sets(list_of_sets: List[set]) -> List[set]:
    """Given list of sets, find sets that intersect and merge them with
    a union operation, then returns the list of sets.  If set does not
    intersect with others, then it is returned without modification.
    """
    output_list = []
    while len(list_of_sets) > 0:
        # Pop sets of indices one by one
        setA = list_of_sets.pop()
        # If first set, add to output list
        if not output_list:
            output_list.append(setA)
        # For later sets, check if set overlaps with existing output_list items
        else:
            intersected = False
            for setB in output_list:
                intersect = setA.intersection(setB)
                # If overlaps, merge them and replace grp in output_list
                if bool(intersect):
                    output = setA.union(setB)
                    output_list.remove(setB)
                    output_list.append(output)
                    intersected = True
            if not intersected:
                output_list.append(setA)
    return output_list


def bootstrapping_mannwhitneyu(data, n_iterations: int = 1000, n_genes: int = 50):
    """
    data= *arrays : sequence of indexables with same length / shape[0]
    Allowed inputs are lists, numpy arrays, scipy-sparse
    matrices or pandas dataframes
    output: lists of pvalues per n_genes for n_iterations
    """
    random_state = 0
    groups = []
    for iteration in range(0, n_iterations):
        res = train_test_split(
            data, random_state=random_state, train_size=0.75, shuffle=True
        )
        groups += [res]

    groups_pvalues = []
    for group in range(0, len(groups)):
        pvalues = []
        for gene in range(0, n_genes):
            __, pvalue = mannwhitneyu(
                groups[group][0][:, gene],
                groups[group][1][:, gene],
                alternative="greater",
            )
            pvalues += [pvalue]
        groups_pvalues += [(pvalues)]
    return groups_pvalues

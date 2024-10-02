from __future__ import print_function, division
from functools import wraps
import numpy as np
from neuron_tree import NeuronTree


class NeuronMorphologyDistribution:
    def __init__(self, distribution):
        assert isinstance(distribution, (np.ndarray, list, tuple))
        if not isinstance(distribution, np.ndarray):
            distribution = np.array(distribution)
        self._distribution = distribution

    def average(self):
        return np.average(self.distribution)

    def min_val(self):
        return np.amin(self.distribution)

    def max_val(self):
        return np.amax(self.distribution)

    def median(self):
        return np.median(self.distribution)

    def variance(self):
        return np.var(self.distribution)

    def std(self):
        return np.std(self.distribution)

    @property
    def distribution(self):
        return self._distribution


class NeuronMorphology:
    def __init__(self, tree):
        assert isinstance(tree, NeuronTree)
        self.tree = tree
        self.tree.find_children(rescan=False)
        self.max_branch_order = max(self.tree.branch_orders().values())
        self.branch_order_tree = None

    def max_branch_order_tree(func):
        def wrapper(self, **kwargs):
            max_branch_order = kwargs.get('max_branch_order')
            max_branch_order = NeuronTree.MAX_BRANCH_ORDER \
                if max_branch_order is None else max_branch_order
            self.branch_order_tree = self.tree \
                if max_branch_order >= self.max_branch_order \
                else self.tree.prune_branch_order(max_branch_order=max_branch_order)
            if len(self.branch_order_tree.children) == 0:
                self.branch_order_tree.find_children()
            func(self, **kwargs)
        return wrapper

    def n_stems(self):
        root_id = self.tree.children[NeuronTree.ROOT_PARENT_ID][0]
        return NeuronMorphologyDistribution(float(len(self.tree.children[root_id])))

    @max_branch_order_tree
    def n_bifurcations(self, max_branch_order=NeuronTree.MAX_BRANCH_ORDER):
        n = 0
        root_id = self.branch_order_tree.children[NeuronTree.ROOT_PARENT_ID][0]
        for node_id, children in self.branch_order_tree.children.items():
            if len(children) > 1 and node_id != root_id:
                n += 1
        return NeuronMorphologyDistribution([float(n)])

    max_branch_order_tree = staticmethod(max_branch_order_tree)


def test():
    swc_path = '/home/muyezhu/mcp3d/python/misc/SW190111-01R_g1_1.swc'
    tree = NeuronTree()
    tree.build_tree_from_swc(swc_path)
    morphology = NeuronMorphology(tree)
    morphology.n_bifurcations(max_branch_order=1)


test()


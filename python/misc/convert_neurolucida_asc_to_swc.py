import os
import argparse
import re
from collections import defaultdict
import json
import numpy as np
from neuron_tree import NeuronTree, NeuronNode


class NeurolucidaTextBlock:
    def __init__(self):
        pass

    @staticmethod
    def n_block_indicators_in_line(line, start=True):
        indicator = '(' if start else ')'
        return line.count(indicator)

    @staticmethod
    def line_block_complete(line):
        return NeurolucidaTextBlock.n_block_indicators_in_line(line, True) > 0 and \
               NeurolucidaTextBlock.n_block_indicators_in_line(line, True) == \
               NeurolucidaTextBlock.n_block_indicators_in_line(line, False)

    @staticmethod
    def line_has_one_complete_block(line):
        return NeurolucidaTextBlock.n_block_indicators_in_line(line, True) == 1 and \
               NeurolucidaTextBlock.n_block_indicators_in_line(line, False) == 1

    @staticmethod
    def next_largest_block(lines):
        pass

    @staticmethod
    def strip_outmost_block_indicators(lines):
        stripped_lines = None
        for i in range(len(lines)):
            if lines[i] == '(':
                if i < len(lines) - 1:
                    stripped_lines = lines[i + 1:]
                    break
                else:
                    return ''
            if i == len(lines) - 1:
                stripped_lines = lines
        for i in range(len(stripped_lines) - 1, -1, -1):
            if stripped_lines[i] == ')':
                stripped_lines = stripped_lines[0:i]
                break
        return stripped_lines


class NeurolucidaNode:
    def __init__(self, x, y, z, r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r

    def __repr__(self):
        output = 'NeurolucidaNode: {} {} {} {}\n'.format(self.x, self.y,
                                                         self.z, self.r)
        return output

    @staticmethod
    def is_node_line(line):
        return NeurolucidaNode.parse_node_line(line) is not None

    @staticmethod
    def parse_node_line(line):
        if not NeurolucidaTextBlock.line_has_one_complete_block(line):
            return None
        line = line.split(')')[0]
        line = line.replace('(', '')
        tokens = line.split()
        if not len(tokens) == 4:
            return None
        try:
            x = float(tokens[0])
            y = float(tokens[1])
            z = float(tokens[2])
            r = float(tokens[3])
        except TypeError:
            return None
        return x, y, z, r


class NeurolucidaContour:
    def __init__(self):
        self.contour_text_block = None
        self.nodes = set()
        self.center = None
        self.radius = None
        self.soma_id = None
        self.node_xyzs = None

    @staticmethod
    def contour_from_text_block(contour_text_block):
        contour = NeurolucidaContour()
        contour.contour_text_block = contour_text_block
        contour.build_contour_from_text()
        return contour

    def __repr__(self):
        output = 'NeurolucidaContour:\n'
        for node in self.nodes:
            output += repr(node)
        return output

    def __iadd__(self, other):
        assert isinstance(other, NeurolucidaContour)
        self.add_nodes(other.nodes)
        if self.node_xyzs is not None and other.node_xyzs is not None:
            self.node_xyzs = np.concatenate((self.node_xyzs, other.node_xyzs),
                                            axis=0)
        return self

    def build_contour_from_text(self):
        soma_pattern = re.compile('Soma ([0-9]+)')
        for l in self.contour_text_block.split('\n'):
            if self.soma_id is None:
                m = re.search(soma_pattern, l)
                if m is not None:
                    self.soma_id = int(m.group(1))
            node_values = NeurolucidaNode.parse_node_line(l)
            if node_values is None:
                continue
            self.add_node(NeurolucidaNode(*node_values))
        assert self.soma_id is not None

    def add_nodes(self, nodes):
        for node in nodes:
            self.add_node(node)

    def add_node(self, node):
        assert isinstance(node, NeurolucidaNode)
        self.nodes.add(node)

    def pool_node_xyzs(self):
        self.node_xyzs = np.zeros(shape=(len(self.nodes), 3))
        for i, node in enumerate(self.nodes):
            self.node_xyzs[i, :] = np.array([node.x, node.y, node.z])

    def calculate_center(self):
        if self.node_xyzs is None:
            self.pool_node_xyzs()
        self.center = np.sum(self.node_xyzs, axis=0) / len(self.nodes)

    def calculate_radius(self):
        if self.node_xyzs is None:
            self.pool_node_xyzs()
        if self.center is None:
            self.calculate_center()
        centers = np.broadcast_to(self.center, shape=self.node_xyzs.shape)
        rs = np.sum(np.square(centers - self.node_xyzs), axis=1)
        self.radius = np.sum(np.sqrt(rs)) / len(self.nodes)

    def distance_to_contour(self, other_node):
        assert isinstance(other_node, NeurolucidaNode)
        if self.node_xyzs is None:
            self.pool_node_xyzs()
        xyz = np.array([other_node.x, other_node.y, other_node.z])
        xyz = np.broadcast_to(xyz, shape=self.node_xyzs.shape)
        rs = np.sum(np.square(xyz - self.node_xyzs), axis=1)
        rs = np.sqrt(rs)
        return np.amin(rs)

    @staticmethod
    def is_contour_block(block):
        return block.find('Soma') >= 0

    @staticmethod
    def merge_contours(contours):
        for contour in contours:
            assert isinstance(contour, NeurolucidaContour)
            if contour.node_xyzs is None:
                contour.pool_node_xyzs()
        merged_contours = {}
        for contour in contours:
            if contour.soma_id not in merged_contours:
                merged_contours[contour.soma_id] = contour
            else:
                merged_contours[contour.soma_id] += contour
        for contour_id, merged_contour in merged_contours.items():
            merged_contour.calculate_center()
            merged_contour.calculate_radius()
        return merged_contours.values()


class NeurolucidaSegment:
    def __init__(self, segment_id, parent_id):
        self.segment_id = segment_id
        self.parent_id = parent_id
        self.nodes = []

    def add_node(self, node):
        assert isinstance(node, NeurolucidaNode)
        self.nodes.append(node)

    def __repr__(self):
        output = 'NeurolucidaSegment: id = {}, parent_id = {}\n'\
                 .format(self.segment_id, self.parent_id)
        for node in self.nodes:
            output += repr(node)
        return output


class NeurolucidaTree:
    def __init__(self, tree_text_block):
        self.tree_text_block = tree_text_block
        self.segment_id = 0
        self.parent_id = -1
        self.segments = defaultdict(lambda: NeurolucidaSegment(self.segment_id,
                                                               self.parent_id))
        self.build_tree()

    def __repr__(self):
        output = 'NeurolucidaTree: {} segments\n'.format(len(self.segments))
        for segment_id, segment in self.segments.items():
            output += repr(segment)
        return output

    def build_tree(self):
        # strip outmost ()
        tree_text_block = NeurolucidaTextBlock.strip_outmost_block_indicators(self.tree_text_block)
        for l in tree_text_block.split('\n'):
            complete_block = NeurolucidaTextBlock.line_block_complete(l)
            node_values = NeurolucidaNode.parse_node_line(l)
            # if contains complete block, but is not a node block, ignore
            if complete_block and node_values is None:
                continue
            # if contains complete block and is a node block
            elif complete_block:
                x, y, z, r = node_values
                self.segments[self.segment_id].add_node(NeurolucidaNode(x, y, z, r))
            # if does not have a complete block
            else:
                # start of children (excluding Varicosity)
                if l.strip() == '(': #
                    self.parent_id = self.segment_id
                    self.segment_id += 1
                # sibling
                elif l.strip() == '|':
                    self.segment_id += 1
                # children end
                elif l.find(')  ;  End of split') >= 0:
                    self.segment_id += 1
                    self.parent_id = self.segments[self.parent_id].parent_id
                # irrelevant characters
                else:
                    continue

    def find_segment_children(self):
        segment_children_ids = defaultdict(list)
        for _, segment in self.segments.items():
            if segment.parent_id == -1:
                continue
            segment_children_ids[segment.parent_id].append(segment.segment_id)
        return segment_children_ids

    def tree_root_node(self):
        return self.segments[0].nodes[0]

    @staticmethod
    def is_tree_block(block):
        return block.find('Dendrite') >= 0


class NeurolucidaASC:
    def __init__(self, asc_path):
        assert os.path.isfile(asc_path)
        self.asc_path = asc_path
        self.contours = []
        self.trees = []

    def build_structures(self):
        for text_block in self.level0_blocks():
            if NeurolucidaContour.is_contour_block(text_block):
                self.contours.append(NeurolucidaContour.contour_from_text_block(text_block))
            elif NeurolucidaTree.is_tree_block(text_block):
                self.trees.append(NeurolucidaTree(text_block))
        self.contours = NeurolucidaContour.merge_contours(self.contours)

    def write_to_swc(self, output_path, connect=True):
        swc_trees = self.extract_neuron_trees(connect=connect)
        out_dir = os.path.dirname(output_path)
        if not os.path.isdir(out_dir) and len(out_dir) > 0:
            os.makedirs(out_dir)
        with open(output_path, 'w') as f:
            for swc_tree_id, swc_tree in swc_trees.items():
                for node in swc_tree.tree.values():
                    NeuronTree.write_swc_file_line(f, node)

    def extract_neuron_trees(self, connect=True):
        swc_trees = defaultdict(NeuronTree)
        swc_tree_id = 0
        node_id = 0
        for contour in self.contours:
            swc_trees[swc_tree_id].add_node(NeuronNode(node_id, 0,
                                                       contour.center[0],
                                                       contour.center[1],
                                                       contour.center[2],
                                                       contour.radius, -1))
            node_id += 1
            swc_tree_id += 1
        for tree in self.trees:
            if connect:
                tree_root_node = tree.tree_root_node()
                distances = [contour.distance_to_contour(tree_root_node)
                             for contour in self.contours]
                swc_tree_id = np.argmin(np.array(distances))
            tree_segment_children = tree.find_segment_children()
            segment_end_node_id = {}
            children = [0]
            while len(children) > 0:
                segment_id = children.pop(-1)
                segment = tree.segments[segment_id]
                children.extend(tree_segment_children[segment_id])
                for i, node in enumerate(segment.nodes):
                    x, y, z, r = node.x, node.y, node.z, node.r
                    if i == 0:
                        if segment_id == 0:
                            if connect:
                                parent_node_id = swc_tree_id
                            else:
                                parent_node_id = NeuronTree.ROOT_PARENT_ID
                        else:
                            parent_node_id = segment_end_node_id[segment.parent_id]
                        # for segments with a single node
                        if i == len(segment.nodes) - 1:
                            segment_end_node_id[segment_id] = node_id
                    elif i == len(segment.nodes) - 1:
                        segment_end_node_id[segment_id] = node_id
                        parent_node_id = node_id - 1
                    else:
                        parent_node_id = node_id - 1
                    swc_trees[swc_tree_id].add_node(NeuronNode(node_id, 0,
                                                               x, y, z, r,
                                                               parent_node_id))
                    node_id += 1
            if not connect:
                swc_tree_id += 1
        return swc_trees

    def write_soma(self, output_json_path):
        output_json = defaultdict(list)
        for contour in self.contours:
            for node in contour.nodes:
                output_json[contour.soma_id].append([node.x, node.y, node.z])
        with open(output_json_path, 'w') as f:
            json.dump(output_json, f, indent=4)

    def print_trees(self):
        for tree in self.trees:
            print(tree)

    def print_contours(self):
        for contour in self.contours:
            print(contour)

    def level0_blocks(self):
        n_block_start_indicators = 0
        n_block_end_indicators = 0
        block_texts = ''
        for l in self._asc_lines():
            block_texts += l
            n_block_start_indicators += NeurolucidaTextBlock.n_block_indicators_in_line(l, True)
            n_block_end_indicators += NeurolucidaTextBlock.n_block_indicators_in_line(l, False)
            if n_block_start_indicators == n_block_end_indicators:
                n_block_start_indicators = 0
                n_block_end_indicators = 0
                texts = block_texts
                block_texts = ''
                yield texts

    def _asc_lines(self):
        with open(self.asc_path, 'r') as f:
            for l in f:
                yield l


def test():
    asc_name = '3_1_bla_lucida.ASC'
    asc = NeurolucidaASC(asc_name)
    asc.build_structures()
    asc.print_trees()
    asc.print_contours()
    asc.write_to_swc('3_1_bla_lucida.swc', False)
    asc.write_soma('3_1_bla_lucida_soma.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--connect', action='store_true',
                        help='connect every neurite to a soma')
    parser.add_argument('input_path', help='path to input swc file')
    parser.add_argument('output_path', help='path to output swc file')
    args = parser.parse_args()
    asc = NeurolucidaASC(args.input_path)
    asc.build_structures()
    asc.write_to_swc(args.output_path, connect=args.connect)
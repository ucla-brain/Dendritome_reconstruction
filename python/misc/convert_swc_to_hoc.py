from __future__ import print_function, division
import os
import sys
sys.path.append('/home/muyezhu/mcp3d/python')
from copy import deepcopy
from collections import defaultdict, namedtuple
from gcut.neuron_tree import NeuronNode
from convert_hoc_to_swc import HocSection, HocFilament

TreeSegment = namedtuple('TreeSegment', ['segment_id', 'segment_node_ids'])


class SwcTreeSegments:
    def __init__(self, swc_trees, swc_node_children):
        self.swc_trees = swc_trees
        self.swc_node_children = swc_node_children
        # segment_start_id: TreeSegment
        self.segments = {}
        self.segment = None
        self.segment_id = None
        self.connections = {}
        self.connection = None
        self.segment_start_ids = None
        self.segment_start_id = None

    def _next_tree_segment(self):
        self.segment_start_id = self.segment_start_ids[0]
        if self.segment_start_id not in self.segments:
            segment = TreeSegment(self.segment_id, [self.segment_start_id])
            self.segments[self.segment_start_id] = segment
            self.segment_id += 1
        else:
            segment = self.segments[self.segment_start_id]

        next_nodes = self.swc_node_children[self.segment_start_id]
        while len(next_nodes) == 1:
            self.segment_start_id = next_nodes[0]
            segment.segment_node_ids.append(self.segment_start_id)
            next_nodes = self.swc_node_children[self.segment_start_id]
        if len(next_nodes) > 1:
            parent_segment_id = self.segments[self.segment_start_ids[0]].segment_id
            self.segment_start_ids.extend(next_nodes)
            for next_node in next_nodes:
                self.segments[next_node] = TreeSegment(self.segment_id, [next_node])
                self.connections[tuple([self.segment_id, 0])] = tuple([parent_segment_id, 1])
                self.segment_id += 1
        self.segment_start_ids.pop(0)

    def tree_segments(self, tree_root_id):
        self.segments = {}
        self.segment_id = None
        self.connections = {}
        self.segment_start_ids = [tree_root_id]
        self.segment_start_id = None
        self.segment_id = 0

        while len(self.segment_start_ids) > 0:
            self._next_tree_segment()

        seen = defaultdict(int)
        for segment in self.segments.values():
            node_ids = segment.segment_node_ids
            for node_id in node_ids:
                seen[node_id] += 1
        for seen_counts in seen.values():
            if seen_counts > 1:
                raise ValueError('node appearing in segments twice')
        for node_id in self.swc_trees[tree_root_id]:
            if node_id not in seen:
                raise ValueError('node not seen in any segments')

        return self.segments, self.connections

    def all_tree_segments(self):
        for tree_root_id in self.swc_trees:
            yield self.tree_segments(tree_root_id)

    def tree_segment_numbers(self, tree_root_id):
        # n_children_root + reduce(number of children if children > 1)
        tree = self.swc_trees[tree_root_id]
        n = 0
        for node_id in tree:
            n_node_children = len(self.swc_node_children[node_id])
            if node_id == tree_root_id or n_node_children > 1:
                n += n_node_children
        return max(n, 1)


class SwcHocConverter:
    Filament_Start_ID = 100000000

    def __init__(self, swc_path):
        if not os.path.isfile(swc_path):
            raise ValueError('{} is not a file'.format(swc_path))
        self.input_path = swc_path
        self.swc_root_ids = set()
        # node_id: NeuronNode
        self.swc_nodes = {}
        # node_id: [children_ids]
        self.swc_node_children = defaultdict(list)
        # root_node_id, {tree_node_ids}
        self.swc_trees = defaultdict(set)
        self.segmentor = None
        self.filaments = []
        self.filament_id = None
        self.output_path = self.input_path.replace('.swc', '.hoc')

    def _feed_line(self):
        with open(self.input_path, 'r') as f:
            for l in f:
                l = l.strip()
                if len(l) == 0 or l[0] == '#':
                    continue
                yield l

    def _filaments_index(self):
        return self.filament_id - SwcHocConverter.Filament_Start_ID

    def extract_swc_nodes(self):
        for line in self._feed_line():
            node_id, node_type, x, y, z, radius, parent_id = line.split()
            if parent_id == '-1':
                self.swc_root_ids.add(int(node_id))
                self.swc_trees[int(node_id)] = set()
            self.swc_nodes[int(node_id)] = \
                NeuronNode(int(node_id), int(node_type),
                           float(x), float(y), float(z),
                           float(radius), int(parent_id))
            if parent_id != '-1':
                self.swc_node_children[int(parent_id)].append(int(node_id))

    def extract_swc_trees(self):
        if len(self.swc_root_ids) == 1:
            swc_root_id = self.swc_root_ids[0]
            self.swc_trees[swc_root_id] = set(self.swc_nodes.keys())
        else:
            for swc_root_id in self.swc_root_ids:
                children_list = deepcopy(self.swc_node_children[swc_root_id])
                while len(children_list) > 0:
                    child = children_list.pop(-1)
                    self.swc_trees[swc_root_id].add(child)
                    children_list.extend(self.swc_node_children[child])
        self.segmentor = SwcTreeSegments(self.swc_trees, self.swc_node_children)

    def build_filaments(self):
        self.filament_id = SwcHocConverter.Filament_Start_ID
        for tree_root_id in self.swc_trees:
            self.add_hoc_filament(tree_root_id)
        for filament in self.filaments:
            assert filament.n_sections == len(filament.sections)

    def write_filament(self, filament_id, f):
        filament = self.filaments[filament_id - SwcHocConverter.Filament_Start_ID]
        f.write(HocFilament.generate_filament_creation_line(
            filament_id, filament.n_sections))
        for section_id in range(filament.n_sections):
            f.write(HocSection.generate_section_begin_line(filament_id, section_id))
            section_nodes = filament.sections[section_id].section_nodes
            for section_node in section_nodes:
                x, y, z, diameter = section_node
                f.write(HocSection.generate_section_node_line(x, y, z, diameter))
            f.write(HocSection.generate_section_end_line())
        for section_connection, parent_section_connection in filament.section_connections.items():
            section_id, section_connector = section_connection
            parent_section_id, parent_section_connector = parent_section_connection
            f.write(HocFilament.generate_filament_section_connection_line(
                filament_id, section_id, section_connector,
                parent_section_id, parent_section_connector))
        f.write(HocFilament.generate_filament_end_line(filament_id))

    def write_filaments(self):
        out_name = os.path.basename(self.output_path).replace('.hoc', '')
        with open(self.output_path, 'w') as f:
            f.write('// converted from {}\n\n'.format(self.input_path))
            f.write('strdef neuron_name\n')
            f.write('neuron_name = \"Filaments 1\"')
            f.write('\n\n')
            for filament_id in range(SwcHocConverter.Filament_Start_ID,
                                     self.filament_id):
                self.write_filament(filament_id, f)

    def add_hoc_section(self, section_id):
        self.filaments[self._filaments_index()].add_section(HocSection(section_id))

    def add_hoc_section_nodes(self, section_id, section_swc_node_ids):
        filament = self.filaments[self._filaments_index()]
        section = filament.sections[section_id]
        for node_id in section_swc_node_ids:
            section.add_node(self.swc_nodes[node_id].x,
                             self.swc_nodes[node_id].y,
                             self.swc_nodes[node_id].z,
                             self.swc_nodes[node_id].radius * 2)

    def add_hoc_section_connection(self, section_id, parent_section_id):
        self.filaments[self._filaments_index()].add_connection(
            self.filament_id, section_id, 0, parent_section_id, 1)

    def add_hoc_filament(self, tree_root_id):
        segments, connections = self.segmentor.tree_segments(tree_root_id)
        self.filaments.append(HocFilament(self.filament_id, len(segments)))
        print(self.filament_id, segments)
        for _, sections in segments.items():
            section_id, section_swc_node_ids = sections
            self.add_hoc_section(section_id)
            self.add_hoc_section_nodes(section_id, section_swc_node_ids)
        for section_connection, parent_connection in connections.items():
            self.add_hoc_section_connection(section_connection[0],
                                            parent_connection[0])
        self.filament_id += 1



swc_path = '1_z257_447_y12881_14073_x5410_6852_cluster_631_run_0_resample_translate.swc'
convertor = SwcHocConverter(swc_path)
convertor.extract_swc_nodes()
convertor.extract_swc_trees()
convertor.build_filaments()
convertor.write_filaments()
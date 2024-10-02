from __future__ import print_function, division
import os
import sys
sys.path.append('/home/muyezhu/mcp3d/python')
from collections import defaultdict
from gcut.neuron_tree import NeuronNode

"""
convert imaris's .hoc neuron morphology files to .swc file
usage: python convert_hoc_to_swc.py hoc_path
output .swc file will be in same directory as input .hoc file and has same 
file name (but with extension replaced)
"""


class HocSection:
    def __init__(self, section_id):
        self.section_id = section_id
        self.section_nodes = []
        self.section_parent_id = None

    def add_node(self, x, y, z, diameter):
        self.section_nodes.append((x, y, z, diameter))

    @staticmethod
    def parse_section_node_line(line):
        # pt3dadd(227.668,395.152,1.46334,3,0)
        if not line.strip().startswith('pt3dadd'):
            return None, None, None, None
        line = line.replace('pt3dadd(', '').replace(')', '')
        x, y, z, diameter, _ = [float(token) for token in line.split(',')]
        return x, y, z, diameter

    @staticmethod
    def parse_filament_section_id(line):
        # filament_100000000[0] {
        if not line.strip().startswith('filament_'):
            return None, None
        line = line.replace('filament_', '').replace('[', ' ').replace(']', '').replace(' {', '')
        filament_id, section_id = [int(token) for token in line.split()]
        return filament_id, section_id

    @staticmethod
    def generate_section_begin_line(filament_id, section_id):
        return 'filament_{}[{}] {{\n  pt3dclear()\n'.format(filament_id, section_id)

    @staticmethod
    def generate_section_node_line(x, y, z, diameter):
        return '  pt3dadd({},{},{},{},0)\n'.format(x, y, z, diameter)

    @staticmethod
    def generate_section_end_line():
        return '}\n\n'


class HocFilament:
    # filament_100000000[1338]
    # 100000000: filament id
    # 1338: section id
    def __init__(self, filament_id, n_sections):
        assert filament_id >= 0 and n_sections > 0
        self.filament_id = filament_id
        self.n_sections = n_sections
        # section_id: HocSection
        self.sections = {}
        # (section_id, 0/1): (section_parent_id, 0/1)
        self.section_connections = {}

    @staticmethod
    def parse_filament_creation(line):
        # create filament_100000000[1339]
        if not line.strip().startswith('create filament'):
            return None, None
        filament_id, n_sections = HocSection.parse_filament_section_id(line.replace('create ', ''))
        return filament_id, n_sections

    @staticmethod
    def parse_section_connection(line):
        if not line.strip().startswith('connect'):
            return None, None, None, None, None
        line = line.replace('connect', '').strip()
        section, parent_section = line.split(', ')
        section_connector = 1 if section.endswith('(1.0)') else 0
        filament_id, section_id = HocSection.parse_filament_section_id(section.replace('(0.0)', '').replace('(1.0)', ''))
        parent_section_connector = 1 if parent_section.endswith('(1.0)') else 0
        _, parent_section_id = HocSection.parse_filament_section_id(parent_section.replace('(0.0)', '').replace('(1.0)', ''))
        assert filament_id == _
        return filament_id, section_id, section_connector, \
               parent_section_id, parent_section_connector

    @staticmethod
    def generate_filament_creation_line(filament_id, n_sections):
        return 'create filament_{}[{}]\n\n'.format(filament_id, n_sections)

    @staticmethod
    def generate_filament_section_connection_line(filament_id, section_id,
                                                  section_connector,
                                                  parent_section_id,
                                                  parent_section_connector):
        return 'connect filament_{}[{}]({}), filament_{}[{}]({})\n'\
                .format(filament_id, section_id,
                        '0.0' if section_connector == 0 else '1.0',
                        filament_id, parent_section_id,
                        '0.0' if parent_section_connector == 0 else '1.0')

    @staticmethod
    def generate_filament_end_line(filament_id):
        return '\ndefine_shape()\naccess filament_{} [0]\n\n'.format(filament_id)

    def add_section(self, section):
        if not isinstance(section, HocSection):
            raise TypeError('expecting HocSection instance')
        if section.section_id in self.sections:
            raise ValueError('attempting to add already existing section')
        self.sections[section.section_id] = section

    def add_connection(self, filament_id, section_id, section_connector,
                       parent_section_id, parent_section_connector):
        assert self.filament_id == filament_id
        assert section_connector == 0 or section_connector == 1
        assert parent_section_connector == 0 or parent_section_connector == 1
        if tuple([section_id, section_connector]) in self.section_connections:
            raise ValueError('section_id {} connecting end {} already has a '
                             'parent'.format(section_id, section_connector))
        self.section_connections[tuple([section_id, section_connector])] = \
            tuple([parent_section_id, parent_section_connector])


class HocToSwcConverter:
    SECTION_NODE_LINE = 0
    SECTION_ID_LINE = 1
    FILAMENT_CREATE_LINE = 2
    SECTION_CONNECT_LINE = 3
    OTHER_LINE = 4

    def __init__(self, input_path):
        if not os.path.isfile(input_path):
            raise ValueError('{} is not a file'.format(input_path))
        self.input_path = input_path
        # filament_id, HocFilament
        self.filaments = {}
        self.active_filament_id = None
        self.active_section_id = None
        self.swc_id = None
        self.swc_nodes = {}
        # (section_id, connector): swc_id
        self.section_connector_swc_ids = {}
        self.output_path = self.input_path.replace('.hoc', '.swc')
        self.actions = {HocToSwcConverter.SECTION_NODE_LINE: self.add_node_to_section,
                        HocToSwcConverter.SECTION_ID_LINE: self.add_section_to_filament,
                        HocToSwcConverter.FILAMENT_CREATE_LINE: self.add_filament,
                        HocToSwcConverter.SECTION_CONNECT_LINE: self.connect_filament_sections,
                        HocToSwcConverter.OTHER_LINE: self._dummy}

    def _feed_line(self):
        with open(self.input_path, 'r') as f:
            for l in f:
                if len(l.strip()) == 0:
                    continue
                yield l

    @staticmethod
    def _geometry_line_type(line):
        line = line.strip()
        if line.startswith('create filament'):
            return HocToSwcConverter.FILAMENT_CREATE_LINE
        elif line.startswith('pt3dadd'):
            return HocToSwcConverter.SECTION_NODE_LINE
        elif line.startswith('filament_') and line.endswith('{'):
            return HocToSwcConverter.SECTION_ID_LINE
        elif line.startswith('connect'):
            return HocToSwcConverter.SECTION_CONNECT_LINE
        else:
            return HocToSwcConverter.OTHER_LINE

    def _dummy(self, line):
        pass

    def convert(self):
        self.build_filaments()
        self.build_swc_trees()
        self.write_swc_trees()

    def build_filaments(self):
        for line in self._feed_line():
            line_type = HocToSwcConverter._geometry_line_type(line)
            if line_type != HocToSwcConverter.OTHER_LINE:
                print(line.strip())
            self.actions[line_type](line)
        for filament_id, filament in self.filaments.items():
            assert filament.n_sections == len(filament.sections)

    def add_filament(self, line):
        filament_id, n_sections = HocFilament.parse_filament_creation(line)
        if filament_id in self.filaments:
            raise ValueError('filament {} already added'.format(filament_id))
        self.filaments[filament_id] = HocFilament(filament_id, n_sections)
        self.active_filament_id = filament_id
        self.active_section_id = None
        print('add filament {}'.format(self.active_filament_id))

    def add_section_to_filament(self, line):
        filament_id, section_id = HocSection.parse_filament_section_id(line)
        assert self.active_filament_id == filament_id
        if filament_id not in self.filaments:
            raise ValueError('attempting to add section to non existing '
                             'filament {}'.format(filament_id))
        self.filaments[filament_id].add_section(HocSection(section_id))
        self.active_section_id = section_id
        print('add section {} to filament {}'.format(self.active_section_id,
                                                     self.active_filament_id))

    def add_node_to_section(self, line):
        x, y, z, diameter = HocSection.parse_section_node_line(line)
        self.filaments[self.active_filament_id].sections[self.active_section_id].add_node(x, y, z, diameter)
        print('add node ({}, {}, {}, {}) to filament {} section {}'
              .format(x, y, z, diameter, self.active_filament_id,
                      self.active_section_id))

    def connect_filament_sections(self, line):
        filament_id, section_id, section_connector, \
        parent_section_id, parent_section_connector = \
            HocFilament.parse_section_connection(line)
        assert self.active_filament_id == filament_id
        self.filaments[self.active_filament_id].add_connection(
            self.active_filament_id, section_id, section_connector,
            parent_section_id, parent_section_connector)
        print('connect filament {} section {} {} to filament {} section {} {}'
              .format(self.active_filament_id, section_id,
                      'start' if section_connector == 0 else 'end',
                      self.active_filament_id, parent_section_id,
                      'start' if parent_section_connector == 0 else 'end'))

    def _section_parent(self, filament_id, section_id):
        filament = self.filaments[filament_id]
        if tuple([section_id, 0]) in filament.section_connections:
            return filament.section_connections[tuple([section_id, 0])]
        elif tuple([section_id, 1]) in filament.section_connections:
            return filament.section_connections[tuple([section_id, 1])]
        else:
            return None

    def _filament_root_section_id(self, filament_id):
        filament = self.filaments[filament_id]
        for section_connect in filament.section_connections:
            section_id, section_connector = section_connect
            break
        while self._section_parent(filament_id, section_id) is not None:
            section_id = self._section_parent(filament_id, section_id)[0]
        return section_id

    def _add_section_to_swc_tree(self, filament, section_id, is_filament_root=False):
        section = filament.sections[section_id]
        for i, node in enumerate(section.section_nodes):
            x, y, z, diameter = node
            radius = diameter / 2
            if 0 < i <= len(section.section_nodes) - 1:
                self.swc_nodes[self.swc_id] = NeuronNode(self.swc_id, 0, x, y, z, radius, self.swc_id - 1)
                if i == len(section.section_nodes) - 1:
                    self.section_connector_swc_ids[tuple([section_id, 1])] = self.swc_id
            else:
                if is_filament_root and i == 0:
                    self.swc_nodes[self.swc_id] = NeuronNode(self.swc_id, 0, x, y, z, radius, -1)
                else:
                    assert tuple([section_id, 0]) in filament.section_connections
                    parent_connection = filament.section_connections[tuple([section_id, 0])]
                    assert parent_connection in self.section_connector_swc_ids
                    self.swc_nodes[self.swc_id] = \
                        NeuronNode(self.swc_id, 0, x, y, z, radius,
                                   self.section_connector_swc_ids[parent_connection])
                self.section_connector_swc_ids[tuple([section_id, 0])] = self.swc_id
            self.swc_id += 1

    def build_swc_tree(self):
        filament = self.filaments[self.active_filament_id]
        # section_id, connector: [(children_id, connector)...]
        reverse_section_connections = defaultdict(list)
        for k, v in filament.section_connections.items():
            reverse_section_connections[v].append(k)
        root_section_id = self._filament_root_section_id(self.active_filament_id)
        self._add_section_to_swc_tree(filament, root_section_id, is_filament_root=True)
        child_sections = reverse_section_connections[tuple([root_section_id, 0])]
        child_sections.extend(reverse_section_connections[tuple([root_section_id, 1])])
        while len(child_sections) > 0:
            section_id, section_connector = child_sections.pop(0)
            self._add_section_to_swc_tree(filament, section_id)
            child_sections.extend(reverse_section_connections[tuple([section_id, 0])])
            child_sections.extend(reverse_section_connections[tuple([section_id, 1])])

    def build_swc_trees(self):
        self.swc_id = 1
        for filament_id in self.filaments:
            self.active_filament_id = filament_id
            self.build_swc_tree()

    def write_swc_trees(self):
        assert len(self.swc_nodes) > 0
        with open(self.output_path, 'w') as f:
            f.write('# converted from {}\n'.format(self.input_path))
            for swc_node in self.swc_nodes.values():
                f.write('{} {} {} {} {} {} {}\n'
                        .format(swc_node.node_id, swc_node.node_type,
                                swc_node.x, swc_node.y, swc_node.z,
                                swc_node.radius, swc_node.parent_id))


def convert(hoc_path):
    converter = HocToSwcConverter(hoc_path)
    converter.convert()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('usage: python convert_hoc_to_swc.py hoc_path')
        exit(1)
    convert(sys.argv[1])

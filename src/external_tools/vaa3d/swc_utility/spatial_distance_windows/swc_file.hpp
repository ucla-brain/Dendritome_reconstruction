//
// Created by muyezhu on 3/3/19.
//

#ifndef MCP3D_VAA3D_SWC_FILE_HPP
#define MCP3D_VAA3D_SWC_FILE_HPP

#include "neuron_tree.hpp"
#include "markers.h"

std::string ReadMetaLines(const std::string &file_path);

// from vaa3d/v3d_external/v3d_main/basic_c_fun/basic_surf_objs.h  mzhu 05/23/2019
NeuronTree readSWC_file_to_tree(const std::string& file_name);

// from vaa3d/v3d_external/v3d_main/basic_c_fun/basic_surf_objs.h  mzhu 05/23/2019
bool writeSWC_file(const std::string& filename, const NeuronTree& nt,
                   const std::string& infostring);

// from vaa3d/vaa3d_tools/released_plugins/v3d_plugins/neurontracing_vn2/app2/my_surf_objs.h  mzhu 05/23/2019
vector<MyMarker*> readSWC_file(string swc_file);
bool readSWC_file(string swc_file, vector<MyMarker> & outmarkers);
bool saveSWC_file(const std::string &swc_path, vector<MyMarker *> &out_markers, const std::string& meta_lines = std::string());
bool saveSWC_file(const std::string &swc_path, vector<MyMarker*> & outmarkers, list<string> & meta_lines);
bool saveSWC_file(const std::string &swc_path, vector<NeuronSWC*> & outmarkers, list<string> & meta_lines);
bool saveSWC_file(const std::string& swc_path, const NeuronTree& neuron_tree, const std::string& meta_lines = std::string());
bool saveDot_file(const std::string &swc_path, vector<MyMarker*> & out_markers);              // save graphviz format

#endif //MCP3D_VAA3D_SWC_FILE_HPP

// Created by muyezhu on 3/2/19.
// this external tool is based on Vaa3D source code by Hanchuan Peng et al
// it reads 2 swc files and computes spatial distance (SD) between 2 neurons
// the swc input file is known to require node index not exceeding total
// node number. it is also assumed each swc file contains a single neuron
// structure. uncertain of other requirements as of yet
// safest to use swc files that have indices from 1 to n, where n is total
// number of nodes in the file, with the indices appearing in ascending order
// looking at the source code, extremely densely sampled swc files will not
// compute properly, resampling is recommended

#include <stdexcept>
#include <fstream>
#include "common/mcp3d_macros.hpp"
#include "vaa3d/swc_utility/neuron_tree.hpp"
#include "vaa3d/swc_utility/swc_file.hpp"
#include "neuron_sim_scores.hpp"

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        MCP3D_MESSAGE("usage: spatial_distance swc_path0 swc_path1 [distance_threshhold] [output_path]")
        return 0;
    }

    double distance_threshold = 2.0;
    if (argc > 3)
    {
        try
        {
            distance_threshold = stof(string(argv[3]));
        }
        catch (const invalid_argument& e)
        {
            MCP3D_MESSAGE("can not convert distance_threshold = " + string(argv[3]) + " to float")
            return 1;
        }
    }
    bool out_file_good = false;
    ofstream out_file;
    if (argc > 4)
    {
        out_file.open(argv[4], ofstream::out);
        out_file_good = out_file.good();
        if (!out_file_good)
        {
            MCP3D_MESSAGE("can not create output file " << argv[4])
            return 1;
        }

    }
    NeuronTree tree0 = readSWC_file_to_tree(string(argv[1]));
    NeuronTree tree1 = readSWC_file_to_tree(string(argv[2]));
    NeuronDistSimple neuron_dist = neuron_score_rounding_nearest_neighbor(&tree0, &tree1, distance_threshold);

    if (!out_file_good)
    {
        cout << "SD_0to1 = " << neuron_dist.dist_12_allnodes << endl;
        cout << "SD_1to0 = " << neuron_dist.dist_21_allnodes << endl;
        cout << "SD = " << neuron_dist.dist_allnodes << endl;
        cout << "SSD = " << neuron_dist.dist_apartnodes << ", distance threshold = " << distance_threshold << endl;
        cout << "SSD% = " << neuron_dist.dist_apartnodes << endl;
        cout << "max distance = " << neuron_dist.dist_max << endl;
    }
    else
    {
        out_file << argv[1] << endl;
        out_file << argv[2] << endl;
        out_file << "SD_0to1 = " << neuron_dist.dist_12_allnodes << endl;
        out_file << "SD_1to0 = " << neuron_dist.dist_21_allnodes << endl;
        out_file << "SD = " << neuron_dist.dist_allnodes << endl;
        out_file << "SSD = " << neuron_dist.dist_apartnodes << ", distance threshold = " << distance_threshold << endl;
        out_file << "SSD% = " << neuron_dist.dist_apartnodes << endl;
        out_file << "max distance = " << neuron_dist.dist_max << endl;
        out_file.close();
    }

    return 0;
}





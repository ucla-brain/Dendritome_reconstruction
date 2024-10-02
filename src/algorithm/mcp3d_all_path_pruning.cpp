//
// Created by muyezhu on 9/20/19.
//
#include <memory>
#include <random>
#include <unordered_map>
#include "image/mcp3d_image_utils.hpp"
#include "mcp3d_vertex.hpp"
#include "mcp3d_all_path_pruning.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    using Vertex = mcp3d::MVertexDijkstra<float>;
    mcp3d::MAllPathPruning<Vertex, uint16_t> app_graph;
    int zdim = 100, ydim = 512, xdim = 512;
    unique_ptr<uint16_t[]> image(new uint16_t[(int64_t)zdim * (int64_t)ydim * (int64_t)xdim]);
    mcp3d::SetRandom<uint16_t>(image.get(), {zdim, ydim, xdim});
    double srcz = 50.0, srcy = 256.0, srcx = 256.0;
    int64_t N = (int64_t)zdim * (int64_t)ydim * (int64_t)xdim;
    //MCP3D_ASSERT(app_graph.ConstructPaths(image.get(), zdim, ydim, xdim, srcz, srcy, srcx, 3))
    app_graph.ConstructPaths(image.get(), zdim, ydim, xdim, srcz, srcy, srcx, 3);
    for (const auto& vertex: app_graph.Graph().Vertices())
    {
        cout << vertex.first << endl;
        MCP3D_ASSERT(vertex.second.has_parent() || vertex.second.is_src())
    }
    app_graph.HierarchicalPrune();
    return 0;
}


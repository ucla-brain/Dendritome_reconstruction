//
// Created by muyezhu on 10/6/19.
//

#include <cmath>
#include <random>
#include <array>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>
#include "common/mcp3d_macros.hpp"
#include "common/mcp3d_utility.hpp"
#include "image/mcp3d_image_utils.hpp"
#include "image/mcp3d_image_in_memory.hpp"
#include "algorithm/mcp3d_vertex.hpp"
#include "algorithm/mcp3d_all_path_pruning.hpp"
#include "test_mcp3d_dijkstra.hpp"

using namespace std;

using AppFloatVertex = mcp3d::MVertexDijkstra<float>;
using AppFloatType = mcp3d::MAllPathPruning<AppFloatVertex, uint16_t>;
using AppFloatGraph = typename AppFloatType::GraphType;
using AppDoubleVertex = mcp3d::MVertexDijkstra<double>;
using AppDoubleType = mcp3d::MAllPathPruning<AppDoubleVertex, uint16_t>;
using AppDoubleGraph = typename AppDoubleType::GraphType;
using TestDoubleVertex = mcp3d::test::TestDijkstraVertex<double>;
using TestDoubleGraph = mcp3d::test::TestDijkstraGraph<TestDoubleVertex>;
using TestDoubleDijkstra = mcp3d::MDijkstra<TestDoubleGraph>;

double VoxelDistanceL2(const vector<int>& dims, int64_t address0, int64_t address1)
{
    array<int, 3> zyx0 = mcp3d::ZyxFromLinearAddress(dims, address0),
                  zyx1 = mcp3d::ZyxFromLinearAddress(dims, address1);
    double deltaz = (double)zyx0[0] - (double)zyx1[0],
           deltay = (double)zyx0[1] - (double)zyx1[1],
           deltax = (double)zyx0[2] - (double)zyx1[2];
    return sqrt(deltaz * deltaz + deltay * deltay + deltax * deltax);
}

template <typename EdgeType, typename VType>
EdgeType VoxelDistanceIntensity(const unique_ptr<VType[]> &data, VType data_max, int64_t vertex_id0, int64_t vertex_id1)
{
    EXPECT_TRUE(data.get());
    VType v0 = data[vertex_id0], v1 = data[vertex_id1];
    EdgeType lambda = (EdgeType)10;
    EdgeType gi0 = exp(lambda *  std::pow(((EdgeType)1.0 -
                                           (EdgeType)v0 / (EdgeType)data_max),
                                          (EdgeType)2)),
             gi1 = exp(lambda *  std::pow(((EdgeType)1.0 -
                                           (EdgeType)v1 / (EdgeType)data_max),
                                          (EdgeType)2));
    return (gi0 + gi1) / 2;
}

double EdgeCost(const unique_ptr<uint16_t[]>& data, uint16_t data_max,
                const vector<int> &dims, int64_t vertex_id0, int64_t vertex_id1)
{
    if (data[vertex_id0] == 0 || data[vertex_id1] == 0)
        return numeric_limits<double>::infinity();
    double d_l2 = VoxelDistanceL2(dims, vertex_id0, vertex_id1),
           d_intensity = VoxelDistanceIntensity<double, uint16_t>(data, data_max, vertex_id0, vertex_id1);
    if (d_l2 > std::sqrt(3.0))
        return std::numeric_limits<double>::infinity();
    return d_l2 * d_intensity;
}

double BernoulliP(const vector<int> &dims, int64_t address, double max_d)
{
    array<int, 3> zyx = mcp3d::ZyxFromLinearAddress(dims, address);
    int64_t center_address = mcp3d::LinearAddress(dims, zyx[0] / 2, zyx[1] / 2, zyx[2] / 2);
    double d = VoxelDistanceL2(dims, address, center_address);
    return 1 - 0.9 * d / max_d;
}

int64_t RandomNeighborAddress(const vector<int> &dims, int64_t address)
{
    array<int, 3> zyx = mcp3d::ZyxFromLinearAddress(dims, address);
    int deltas[3] = {-1, 0, 1};
    int64_t seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<int> distribution(0, 2);
    int delta;
    int zyx_neighbor[3];
    for (int i = 0; i < 3; ++i)
    {
        do
        {
            delta = deltas[distribution(generator)];
            zyx_neighbor[i] = zyx[i] + delta;
        } while (zyx_neighbor[i] < 0 || zyx_neighbor[i] >= dims[i]);
    }
    return mcp3d::LinearAddress(dims, zyx_neighbor[0], zyx_neighbor[1], zyx_neighbor[2]);
}

int64_t RandomAddress(const vector<int64_t>& addresses)
{
    int64_t seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<size_t> uniform(0, addresses.size() - 1);
    size_t target_position = uniform(generator);
    return *(addresses.begin() + target_position);
}

template <typename VType>
void GrowSparseVolume(std::unique_ptr<VType []>& data, double percent,
                      int zdim, int ydim, int xdim)
{
    std::vector<int> dims = {zdim, ydim, xdim};
    // allocate volume and set all value to zero
    size_t n_total_voxels = (size_t)zdim * (size_t)ydim * (size_t)xdim;
    data = std::make_unique<VType []>(n_total_voxels);
    auto data_ptr = data.get();
    mcp3d::SetConstant(data_ptr, dims, (uint16_t)0);
    // set foreground voxel number and source vertex
    auto n_fg_voxels = (int64_t)(round((double)n_total_voxels * percent));
    int srcz = zdim / 2, srcy = ydim / 2, srcx = xdim / 2;
    int64_t address = mcp3d::LinearAddress(dims, srcz, srcy, srcx);
    data_ptr[address] = 1;
    int64_t neighbor_address = RandomNeighborAddress(dims, address);
    vector<int64_t> fg_addresses({address});

    // generate random foreground voxels until number of foreground voxels
    // is equal to n_fg_voxels
    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<uint16_t> uniform(1, numeric_limits<uint16_t>::max());

    // bernoulli random variable: p decreases moving away from the source
    double max_d = 0.5 * VoxelDistanceL2(dims, 0, mcp3d::ReduceProdSeq<int64_t>(dims) - 1);
    while (fg_addresses.size() < (size_t)n_fg_voxels)
    {
        // identify a background voxel neighboring an existing foreground voxel
        while (data_ptr[neighbor_address] > 0)
        {
            address = RandomAddress(fg_addresses);
            neighbor_address = RandomNeighborAddress(dims, address);
        }
        // dice roll on if to set the background voxel to foreground
        double p = BernoulliP(dims, neighbor_address, max_d);
        bernoulli_distribution bernoulli(p);
        bool grow = bernoulli(generator);
        if (grow)
        {
            data_ptr[neighbor_address] = uniform(generator);
            fg_addresses.push_back(neighbor_address);
        }
        address = RandomAddress(fg_addresses);
        neighbor_address = RandomNeighborAddress(dims, address);
    }
}

void GenerateTestDijkstraGraph(TestDoubleGraph& test_graph,
                               unique_ptr<uint16_t[]>& data, const vector<int> &dims)
{
    auto data_ptr = data.get();
    auto data_max = mcp3d::VolumeMax<uint16_t>(data.get(), dims);
    // create vertices
    auto n_total = mcp3d::ReduceProdSeq<int64_t>(dims);
    for (int64_t i = 0; i < n_total; ++i)
        if (data_ptr[i] > 0)
            test_graph[i] = TestDoubleVertex(i);
    // create edges. compute edge cost same as mcp3d::MGraphAllPath
    unique_ptr<int64_t []> candidate_neighbors;
    int64_t n_neighbors;
    for (auto& item: test_graph.Vertices())
    {
        int64_t vertex_id = item.first;
        TestDoubleVertex& vertex = item.second;
        mcp3d::NeighborAddresses(dims, vertex_id, 3, candidate_neighbors);
        n_neighbors = candidate_neighbors[0];
        for (int64_t i = 1; i < n_neighbors; ++i)
        {
            int64_t neighbor_vertex_id = candidate_neighbors[i];
            if (data_ptr[neighbor_vertex_id] == 0)
                continue;
            EXPECT_TRUE(test_graph.HasVertex(neighbor_vertex_id));
            unordered_map<int64_t, double>& neighbor_edges = test_graph.At(neighbor_vertex_id).edges_ref();
            if (neighbor_edges.find(vertex_id) != neighbor_edges.end())
                vertex.edges_ref()[neighbor_vertex_id] = neighbor_edges[vertex_id];
            else
            {
                double edge_cost = EdgeCost(data, data_max, dims, vertex_id, neighbor_vertex_id);
                vertex.edges_ref()[neighbor_vertex_id] = edge_cost;
                neighbor_edges[vertex_id] = edge_cost;
            }
        }
    }
}

// compare against equivalent problem solved using test_dijkstra graph
// edge cost is of double type. source node is always at center of volume
TEST(MAllPathPruning, ConstructPaths)
{
    AppDoubleType app_double;
    AppDoubleGraph& app_graph = app_double.Graph();
    TestDoubleGraph test_graph;
    TestDoubleDijkstra test_dijkstra;

    // create a random sparse problem in 3d array, translate the problem
    // to an equivalent TestDijkstraGraph
    double percent = 0.03;
    int zdim = 256, ydim = 256, xdim = 256;
    int srcz = zdim / 2, srcy = ydim / 2, srcx = xdim / 2;
    int64_t source_id = mcp3d::LinearAddress({zdim, ydim, xdim}, srcz, srcy, srcx);
    unique_ptr<uint16_t[]> data;
    GrowSparseVolume<uint16_t>(data, percent, zdim, ydim, xdim);
    app_graph.Init(data.get(), zdim, ydim, xdim, srcz, srcy, srcx, 3);
    GenerateTestDijkstraGraph(test_graph, data, {zdim, ydim, xdim});
    // validate equivalence of represented problem in app_graph and test_graph
    // equal number of vertices
    EXPECT_EQ(test_graph.Size(), app_graph.Size());
    // same set of vertex ids
    for (const auto& item: app_graph.Vertices())
        EXPECT_TRUE(test_graph.Vertices().find(item.first) != test_graph.Vertices().cend());
    // same edges for each vertex
    vector<int64_t> app_neighbor_ids, test_neighbor_ids;
    for (const auto& item: app_graph.Vertices())
    {
        int64_t vertex_id = item.first;
        EXPECT_TRUE(app_graph.HasVertex(vertex_id));
        app_neighbor_ids = app_graph.NeighborVertexIds(vertex_id);
        test_neighbor_ids = test_graph.NeighborVertexIds(vertex_id);
        sort(app_neighbor_ids.begin(), app_neighbor_ids.end());
        sort(test_neighbor_ids.begin(), test_neighbor_ids.end());
        EXPECT_EQ(test_neighbor_ids, app_neighbor_ids);
        for (const auto& neighbor_id: app_neighbor_ids)
            EXPECT_DOUBLE_EQ(test_graph.EdgeCost(vertex_id, neighbor_id), app_graph.EdgeCost(vertex_id, neighbor_id));
    }
    // dijkstra's shortest path on test_graph and app_graph
    INIT_TIMER(0)
    TIC(0)
    app_double.ConstructPaths(data.get(), zdim, ydim, xdim, srcz, srcy, srcx, 3);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("app", 0)
    TIC(0)
    test_dijkstra.FindPaths(test_graph, source_id);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("dijkstra", 0)
    for (const auto& item: app_graph.Vertices())
    {
        int64_t vertex_id = item.first;
        // the solution should for each vertex have either identical parent
        // or if parent is different, their distance to source node are the same
        if (test_graph.ParentVertexId(vertex_id) != app_graph.ParentVertexId(vertex_id))
            EXPECT_DOUBLE_EQ(test_graph.Dsrc(vertex_id), app_graph.Dsrc(vertex_id));
    }
}


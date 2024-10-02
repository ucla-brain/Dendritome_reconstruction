//
// Created by muyezhu on 10/15/19.
//
#include <type_traits>
#include <limits>
#include <memory>
#include <chrono>
#include <cmath>
#include <random>
#include <gtest/gtest.h>
#include "common/mcp3d_exceptions.hpp"
#include "common/mcp3d_utility.hpp"
#include "image/mcp3d_image_utils.hpp"
#include "algorithm/mcp3d_vertex.hpp"
#include "algorithm/mcp3d_graph.hpp"

using namespace std;

TEST(MGraphAllPath, Constructor)
{
    using Vertex = mcp3d::MVertexDijkstra<double>;
    using GraphDoubleData = mcp3d::MGraphAllPath<Vertex, double>;
    GraphDoubleData graph_double;
    bool correct_vertex_type = is_same<Vertex, GraphDoubleData::VertexType>();
    EXPECT_TRUE(correct_vertex_type);
    bool correct_edge_type = is_same<double, GraphDoubleData::EdgeType>();
    EXPECT_TRUE(correct_edge_type);
    bool correct_data_type = is_same<double, GraphDoubleData::DataType>();
    EXPECT_TRUE(correct_data_type);
    EXPECT_EQ(vector<int>({0, 0, 0}), graph_double.data_dims());
    EXPECT_EQ(numeric_limits<int64_t>::lowest(), graph_double.SrcId());
    EXPECT_EQ(3, graph_double.cnn_type());
    EXPECT_EQ(0, graph_double.n_total_voxels());
    // thrown by mcp3d::MGraphAllPath<Vertex, double>::At
    EXPECT_THROW(graph_double.VertexValue(0), mcp3d::MCPOutOfRangeError);
    EXPECT_EQ((size_t)0, graph_double.Size());
}

// create data volume. set half of the volume to zero:
// even number indices are zero, odd number indices are positive
void PopulateVolume(unique_ptr<double []>& data, int zdim, int ydim, int xdim)
{
    auto n = mcp3d::ReduceProdSeq<int64_t>(vector<int>({zdim, ydim, xdim}));
    ASSERT_TRUE(n > 0);
    data = make_unique<double []>((size_t)n);
    double* data_ptr = data.get();

    mcp3d::SetRandom<double>(data.get(), {zdim, ydim, xdim});
    for (int64_t i = 0; i < n; ++i)
    {
        if (i % 2 == 0)
            data_ptr[i] = (double)0;
        else if (abs(data_ptr[i]) <= 1e-6)
            data_ptr[i] += 1.0;
        else if (data_ptr[i] < 0)
            data_ptr[i] = - data_ptr[i];
    }
}

TEST(MGraphAllPath, Init)
{
    using Vertex = mcp3d::MVertexDijkstra<double>;
    using GraphDoubleData = mcp3d::MGraphAllPath<Vertex, double>;
    GraphDoubleData graph_double;

    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    int zdim = 512, ydim = 512, xdim = 512;
    auto n = mcp3d::ReduceProdSeq<int64_t>(vector<int>({zdim, ydim, xdim}));
    unique_ptr<double []> data;
    double* data_ptr = data.get();
    // note that source vertex is even
    int srcz = zdim / 2, srcy = ydim / 2, srcx = xdim / 2;
    int64_t source_id = mcp3d::LinearAddress({zdim, ydim, xdim}, srcz, srcy, srcx);
    // null data ptr
    EXPECT_THROW(graph_double.Init(data_ptr, zdim, ydim, xdim, srcz, srcy, srcx, 1), mcp3d::MCPRuntimeError);
    PopulateVolume(data, zdim, ydim, xdim);
    // non positve source voxel
    data_ptr = data.get();
    EXPECT_THROW(graph_double.Init(data_ptr, zdim, ydim, xdim, srcz, srcy, srcx, 1), mcp3d::MCPRuntimeError);
    // set source voxel to positive
    data_ptr[source_id] += 1.0;
    graph_double.Init(data_ptr, zdim, ydim, xdim, srcz, srcy, srcx, 1);
    // number of vertices: number of odd number indices + 1 (source id is positive while even)
    EXPECT_EQ((size_t)(n / 2) + 1, graph_double.Size());
    // vertex ids correct
    for (const auto& item: graph_double.Vertices())
        EXPECT_TRUE(item.first % 2 > 0 || item.first == source_id);
    // data dimension correct
    EXPECT_EQ(vector<int>({zdim, ydim, xdim}), graph_double.data_dims());
    // total voxels correct
    EXPECT_EQ(n, graph_double.n_total_voxels());
    // source id correct
    EXPECT_EQ(source_id, graph_double.SrcId());
    // vertex values and dsrc values correct
    auto data_max = mcp3d::VolumeMax<double>(data_ptr, {zdim, ydim, xdim});
    uniform_int_distribution<int64_t> distribution(0, n - 1);
    for (int i = 0; i < 1000; ++i)
    {
        int64_t address = distribution(generator);
        if (address % 2 == 0 && address != source_id)
            EXPECT_FALSE(graph_double.HasVertex(address));
        else
        {
            double expeted_value = exp(10.0 *  std::pow((1.0 - data_ptr[address] / data_max), 2.0));
            EXPECT_DOUBLE_EQ(expeted_value, graph_double.VertexValue(address));
            EXPECT_EQ(numeric_limits<double>::infinity(), graph_double.Dsrc((address)));
        }
    }
    // cnn type correct
    EXPECT_EQ(1, graph_double.cnn_type());

    // re initialize with new data
    zdim = ydim = xdim = 256;
    // note that source vertex is odd and positive
    srcz = zdim / 2;
    srcy = ydim / 2;
    srcx = xdim / 2 + 1;
    source_id = mcp3d::LinearAddress({zdim, ydim, xdim}, srcz, srcy, srcx);
    n = mcp3d::ReduceProdSeq<int64_t>(vector<int>({zdim, ydim, xdim}));
    PopulateVolume(data, zdim, ydim, xdim);
    data_ptr = data.get();
    // default cnn type = 3
    graph_double.Init(data_ptr, zdim, ydim, xdim, srcz, srcy, srcx);
    EXPECT_EQ((size_t)(n / 2), graph_double.Size());
    // vertex ids correct
    for (const auto& item: graph_double.Vertices())
        EXPECT_TRUE(item.first % 2 > 0 || item.first == source_id);
    // data dimension correct
    EXPECT_EQ(vector<int>({zdim, ydim, xdim}), graph_double.data_dims());
    // total voxels correct
    EXPECT_EQ(n, graph_double.n_total_voxels());
    // source id correct
    EXPECT_EQ(source_id, graph_double.SrcId());
    // vertex values and dsrc values correct
    data_max = mcp3d::VolumeMax<double>(data_ptr, {zdim, ydim, xdim});
    distribution = uniform_int_distribution<int64_t>(0, n - 1);
    for (int i = 0; i < 1000; ++i)
    {
        int64_t address = distribution(generator);
        if (address % 2 == 0 && address != source_id)
            EXPECT_FALSE(graph_double.HasVertex(address));
        else
        {
            double expeted_value = exp(10.0 *  std::pow((1.0 - data_ptr[address] / data_max), 2.0));
            EXPECT_DOUBLE_EQ(expeted_value, graph_double.VertexValue(address));
            EXPECT_EQ(numeric_limits<double>::infinity(), graph_double.Dsrc((address)));
        }
    }
    // cnn type correct
    EXPECT_EQ(3, graph_double.cnn_type());
}

TEST(MGraphAllPath, ParentVertex)
{
    using Vertex = mcp3d::MVertexDijkstra<double>;
    using GraphDoubleData = mcp3d::MGraphAllPath<Vertex, double>;
    GraphDoubleData graph_double;

    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    int zdim = 512, ydim = 512, xdim = 512;
    auto n = mcp3d::ReduceProdSeq<int64_t>(vector<int>({zdim, ydim, xdim}));
    unique_ptr<double []> data;
    double* data_ptr = data.get();
    // note that source vertex is odd and positive
    int srcz = zdim / 2, srcy = ydim / 2, srcx = xdim / 2 + 1;
    PopulateVolume(data, zdim, ydim, xdim);
    // non positve source voxel
    data_ptr = data.get();
    graph_double.Init(data_ptr, zdim, ydim, xdim, srcz, srcy, srcx, 3);
    // no vertex has parent
    uniform_int_distribution<int64_t> distribution(0, n - 1);
    for (const auto& item: graph_double.Vertices())
    {
        EXPECT_EQ(numeric_limits<int64_t>::lowest(), graph_double.ParentVertexId(item.first));
        EXPECT_EQ(nullptr, graph_double.ParentVertex(item.first));
    }
    // setting parent vertex correct
    unique_ptr<int64_t []> neighbor_addresses;
    for (int i = 0; i < 100; ++i)
    {
        int64_t address = 0, parent_address = 0;
        do
        {
            address = distribution(generator);
        } while (abs(data_ptr[address]) <= 1e-6);
        mcp3d::NeighborAddresses({zdim, ydim, xdim}, address, 3, neighbor_addresses);
        for (int j = 1; j < neighbor_addresses[0]; ++j)
        {
            parent_address = neighbor_addresses[j];
            if (abs(data_ptr[address]) > 1e-6)
                break;
        }
        graph_double.SetParentVertexId(address, parent_address);
        EXPECT_EQ(parent_address, graph_double.ParentVertexId(address));
    }
}

TEST(MGraphAllPath, NeighborVertexIds)
{
    using Vertex = mcp3d::MVertexDijkstra<double>;
    using GraphDoubleData = mcp3d::MGraphAllPath<Vertex, double>;
    GraphDoubleData graph_double;

    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    int zdim = 512, ydim = 512, xdim = 512 + 1;
    auto n = mcp3d::ReduceProdSeq<int64_t>(vector<int>({zdim, ydim, xdim}));
    uniform_int_distribution<int64_t> distribution(0, n - 1);
    unique_ptr<double []> data;
    double* data_ptr = data.get();
    // note that source vertex is odd and is populated to zero
    int srcz = zdim / 2, srcy = ydim / 2, srcx = xdim / 2 + 1;

    PopulateVolume(data, zdim, ydim, xdim);
    data_ptr = data.get();
    // connectivity_type = 1
    int connectivity_type = 1;
    graph_double.Init(data_ptr, zdim, ydim, xdim, srcz, srcy, srcx, connectivity_type);
    int64_t address;
    vector<int64_t> graph_neighbors;
    unique_ptr<int64_t []> candidate_neighbors;
    unordered_set<int64_t> expected_neighbors;
    for (int i = 0; i < 100; ++i)
    {
        expected_neighbors.clear();
        do
        {
            address = distribution(generator);
        } while (abs(data_ptr[address]) <= 1e-6);
        graph_neighbors = graph_double.NeighborVertexIds(address);
        mcp3d::NeighborAddresses({zdim, ydim, xdim}, address, connectivity_type, candidate_neighbors);
        for (int j = 1; j <= candidate_neighbors[0]; ++j)
            if (abs(data_ptr[candidate_neighbors[j]]) > 1e-6)
                expected_neighbors.insert(candidate_neighbors[j]);
        EXPECT_EQ(expected_neighbors.size(), graph_neighbors.size());
        for (const auto& graph_neighbor: graph_neighbors)
            EXPECT_TRUE(expected_neighbors.find(graph_neighbor) != expected_neighbors.end());
    }
    // connectivity_type = 3
    connectivity_type = 3;
    // all voxels are positve. therefore all spatially valid neighbors are valid
    // neighbors
    mcp3d::SetConstant<double>(data_ptr, {zdim, ydim, xdim}, 1.0);
    graph_double.Init(data_ptr, zdim, ydim, xdim, srcz, srcy, srcx, connectivity_type);
    for (int i = 0; i < 100; ++i)
    {
        do
        {
            address = distribution(generator);
        } while (abs(data_ptr[address]) <= 1e-6);
        graph_neighbors = graph_double.NeighborVertexIds(address);
        mcp3d::NeighborAddresses({zdim, ydim, xdim}, address, connectivity_type, candidate_neighbors);
        EXPECT_EQ((size_t)candidate_neighbors[0], graph_neighbors.size());
    }
}


//
// Created by muyezhu on 10/1/19.
//
#include <cstdint>
#include <string>
#include <fstream>
#include <type_traits>
#include <limits>
#include <utility>
#include <chrono>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <gtest/gtest.h>
#include "common/mcp3d_paths.hpp"
#include "common/mcp3d_utility.hpp"
#include "algorithm/mcp3d_vertex.hpp"
#include "algorithm/mcp3d_graph.hpp"
#include "algorithm/mcp3d_dijkstra.hpp"
#include "test_mcp3d_dijkstra.hpp"

using namespace std;

using TestDoubleVertex = mcp3d::test::TestDijkstraVertex<double>;
using TestDoubleGraph = mcp3d::test::TestDijkstraGraph<TestDoubleVertex>;
using VertexPath = mcp3d::test::VertexPath;

// correctness of TestDijkstraGraph<Vertex>::AllPaths using known static graph
TEST(TestDijkstraGraph, FromFileAllPaths)
{
    string vertex_file_path0(mcp3d::JoinPath(mcp3d::test_data_dir(), "algorithm", "TestDijkstraGraph_fromfile0"));
    TestDoubleGraph graph0(vertex_file_path0);
    EXPECT_EQ((size_t)8, graph0.Size());
    vector<VertexPath> all_paths = graph0.AllPaths(0, 7);
    vector<size_t> path_lengths;
    for (const auto& path: all_paths)
        path_lengths.push_back(path.NVertices());
    sort(path_lengths.begin(), path_lengths.end());
    EXPECT_EQ(vector<size_t>({4, 5, 5, 5, 6, 6}), path_lengths);
    // no paths from vertex 3 to vertex 0
    EXPECT_TRUE(graph0.AllPaths(3, 0).empty());

    path_lengths.clear();
    string vertex_file_path1(mcp3d::JoinPath(mcp3d::test_data_dir(), "algorithm", "TestDijkstraGraph_fromfile1"));
    TestDoubleGraph graph1(vertex_file_path1);
    EXPECT_EQ((size_t)6, graph1.Size());
    all_paths = graph1.AllPaths(0, 4);
    for (const auto& path: all_paths)
        path_lengths.push_back(path.NVertices());
    sort(path_lengths.begin(), path_lengths.end());
    EXPECT_EQ(vector<size_t>({3, 3, 4, 4, 4, 4, 5, 5, 5, 6}), path_lengths);
    // no paths from vertex 5 to vertex 3
    EXPECT_TRUE(graph0.AllPaths(5, 3).empty());
}

// n * (n - 1) * ... * (n - m + 1)
int64_t Permutation(int n, int m)
{
    EXPECT_TRUE(n >= m);
    int64_t result = 1;
    for (int64_t i = 0; i < (int64_t)m; ++i)
        result *= ((int64_t)n - i);
    return result;
}

int64_t PermutationSum(int n)
{
    EXPECT_TRUE(n >= 0);
    int64_t result = 0;
    for (int m = 0; m <= n; ++m)
        result += Permutation(n, m);
    return result;
}

// correctness of TestDijkstraGraph<Vertex>::AllPaths using random graph
// the random graphs are fully connected, therefore many properties are known
TEST(TestDijkstraGraph, RandomGraphAllPaths)
{
    TestDoubleGraph random_graph;
    vector<VertexPath> all_paths;
    vector<size_t> path_edges;
    vector<int> n_vertices({1, 2, 3, 4, 5, 6, 7, 8});
    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    for (int n_vertice: n_vertices)
    {
        // generate fully connected graphs: each vertex has an out edge to
        // all vertices in graph (for graph with vertices 0, 1, 2, 3, 4,
        // vertex 4 has 5 out edges 40, 41, 42, 43, 44)
        // in such graphs, the total number of paths between any two vertices
        // are known, as well as the number of edges in paths with minimally
        // and maximally possible number of edges
        int n_edges = n_vertice * n_vertice;
        cout << "fully connected graph: " << n_vertice << " vertices, " << n_edges << " edges" << endl;
        random_graph.GenerateRandomGraph(n_vertice, n_edges, false);
        EXPECT_EQ((size_t)n_vertice, random_graph.Size());
        EXPECT_EQ((size_t)n_edges, random_graph.NEdges());
        // take random source and target id
        uniform_int_distribution<int64_t> vertex_distribution((int64_t)0, (int64_t)(n_vertice - 1));
        int64_t source_id = vertex_distribution(generator),
                target_id = vertex_distribution(generator);;

        // find all paths between source and target
        all_paths = random_graph.AllPaths(source_id, target_id);
        // if source_id != target_id, [0, n_vertice - 2] number of vertices
        // can be inserted between source and target. if source_id == target_id,
        // number of vertices in between lies in [0, n_vertice - 1]
        EXPECT_EQ(PermutationSum(n_vertice - (source_id == target_id ? 1 : 2)),
                  (int64_t)all_paths.size());

        for (const auto& path: all_paths)
            path_edges.push_back(path.NEdges());
        sort(path_edges.begin(), path_edges.end());

        EXPECT_EQ((size_t)1, path_edges[0]);
        EXPECT_EQ((size_t)(source_id == target_id ? n_vertice : n_vertice - 1),
                  path_edges[path_edges.size() - 1]);
        path_edges.clear();
    }
}

// correctness of MDijkstra<Graph>'s aliases
TEST(MDijkstra, GraphType)
{
    mcp3d::MDijkstra<TestDoubleGraph> dijkstra;
    using DijkstraType = mcp3d::MDijkstra<TestDoubleGraph>;
    bool same_graph_type = is_same<typename DijkstraType::GraphType, TestDoubleGraph>();
    EXPECT_TRUE(same_graph_type);
}

// correctness of MDijkstra<Graph>'s constructor
TEST(MDijkstra, Constructor)
{
    mcp3d::MDijkstra<TestDoubleGraph> dijkstra;
    EXPECT_TRUE(dijkstra.target_ids().empty());
    EXPECT_TRUE(dijkstra.unreached_ids().empty());
    EXPECT_EQ(numeric_limits<int64_t>::lowest(), dijkstra.source_id());
    EXPECT_FALSE(dijkstra.track_boundary());
    EXPECT_FALSE(dijkstra.find_path_complete());
}

template <typename Graph>
pair<double, vector<int64_t>> ShortestOfAllPaths(const Graph& graph, const vector<VertexPath>& all_paths)
{
    if (all_paths.empty())
        return make_pair(numeric_limits<double>::infinity(), vector<int64_t>{});
    vector<double> dsrcs;
    for (const auto& path: all_paths)
    {
        double dsrc = 0;
        for (size_t i = 0; i < path.NEdges(); ++i)
            dsrc += graph.EdgeCost(path.vertex_sequence()[i],
                                   path.vertex_sequence()[i + 1]);
        dsrcs.push_back(dsrc);
    }
    int64_t shortest_path_index = min_element(dsrcs.begin(), dsrcs.end()) - dsrcs.begin();
    VertexPath shortest_path = all_paths[shortest_path_index];
    vector<int64_t> shortest_path_vec(shortest_path.vertex_sequence());
    return make_pair(dsrcs[shortest_path_index], shortest_path_vec);
}

// correctness of MDijkstra<Graph> path finding from source to single target
// the generated path. all possible paths between source and target are found
// by TestDijkstraGraph::AllPath. the path with smallest cost is compared with
// dijkstra's result
TEST(MDijkstra, FindPathSingleTarget)
{
    using DijkstraType = mcp3d::MDijkstra<TestDoubleGraph>;
    DijkstraType dijkstra;
    typename DijkstraType::GraphType test_graph;
    vector<int> n_vertices({10, 20});
    vector<double> n_edge_percents({0.025, 0.05, 0.1, 0.2});
    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    for (int n_vertice: n_vertices)
    {
        // generate random graph with specified number of vertices and edges
        for (double n_edge_percent: n_edge_percents)
        {
            auto n_edge = (int)(n_edge_percent * n_vertice * n_vertice);
            cout << "random test case: " << n_vertice << " vertices, " << n_edge << " edges" << endl;
            test_graph.GenerateRandomGraph(n_vertice, n_edge);
            EXPECT_EQ((size_t)n_vertice, test_graph.Size());
            EXPECT_EQ((size_t)n_edge, test_graph.NEdges());
            // take random source and target id
            uniform_int_distribution<int64_t> vertex_distribution((int64_t)0, (int64_t)(n_vertice - 1));
            int64_t source_id, target_id;
            do
            {
                source_id = vertex_distribution(generator);
                target_id = vertex_distribution(generator);
            } while (source_id == target_id);
            // find all paths between source and target
            cout << "finding all paths from source and target" << endl;
            vector<VertexPath> all_paths = test_graph.AllPaths(source_id, target_id);
            // shortest_path_info: smallest dsrc, cooresponding vector of vertex id
            auto shortest_path_info = ShortestOfAllPaths<typename DijkstraType::GraphType>(test_graph, all_paths);
            // if ant paths exist between and source and target, the shortest such path should contain
            // at least two vertices
            EXPECT_TRUE(shortest_path_info.second.size() >= (size_t)2 || shortest_path_info.second.empty());
            cout << "finding shortest paths from source vertex id = " << source_id << " to target vertex id = " << target_id << endl;
            bool has_path = dijkstra.FindPaths(test_graph, source_id, {target_id}, false);
            EXPECT_TRUE(dijkstra.find_path_complete());
            EXPECT_EQ(unordered_set<int64_t>({target_id}), dijkstra.target_ids());
            EXPECT_TRUE(test_graph.ParentVertexId(source_id) < 0);
            if (has_path)
            {
                cout << "shortest path found" << endl;
                EXPECT_EQ(shortest_path_info.first, test_graph.Dsrc(target_id));
                cout << "shortest distance = " << test_graph.Dsrc(target_id) << endl;
                int64_t current_id = target_id;
                for (auto rit = shortest_path_info.second.crbegin();
                     rit != shortest_path_info.second.crend(); ++rit)
                {
                    EXPECT_EQ(*rit, current_id);
                    current_id = test_graph.ParentVertexId(current_id);
                }
                cout << "shortest path = " << mcp3d::JoinVector(shortest_path_info.second) << endl;
            }
            else
            {
                EXPECT_TRUE(all_paths.empty());
                cout << "no path found" << endl;
            }
        }
    }
}

void MDijkstraFindPathMultiTargetsImpl()
{
    using DijkstraType = mcp3d::MDijkstra<TestDoubleGraph>;
    DijkstraType dijkstra;
    typename DijkstraType::GraphType test_graph;
    int n_vertices = 20;
    double n_edge_percent = 0.15;
    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    auto n_edge = (int)(n_edge_percent * n_vertices * n_vertices);
    cout << "random test case: " << n_vertices << " vertices, " << n_edge << " edges" << endl;
    test_graph.GenerateRandomGraph(n_vertices, n_edge);
    // take random source id
    uniform_int_distribution<int64_t> vertex_distribution((int64_t)0, (int64_t)(n_vertices - 1));
    int64_t source_id = vertex_distribution(generator);
    // find shortest paths between all target and source
    dijkstra.FindPaths(test_graph, source_id, {}, false);
    unordered_set<int64_t> target_ids;
    for (int64_t target_id = 0; target_id < n_vertices; ++target_id)
        if (target_id != source_id)
            target_ids.insert(target_id);
    EXPECT_EQ(target_ids, dijkstra.target_ids());
    EXPECT_TRUE(test_graph.ParentVertexId(source_id) < 0);
    // for each target_id, find all paths between source and target
    unordered_map<int64_t, vector<VertexPath>> target_all_paths;
    for (int64_t target_id = 0; target_id < n_vertices; ++target_id)
    {
        if (target_id == source_id)
            continue;
        target_all_paths[target_id] = test_graph.AllPaths(source_id, target_id);
    }
    unordered_set<int64_t> reached_ids, unreached_ids;
    for (int64_t target_id = 0; target_id < n_vertices; ++target_id)
    {
        if (target_id == source_id)
            continue;
        if (target_all_paths[target_id].empty())
            unreached_ids.insert(target_id);
        else
            reached_ids.insert(target_id);
    }
    // correctness of reached and unreached targets
    EXPECT_EQ(reached_ids, dijkstra.ReachedTargetIds());
    EXPECT_EQ(unreached_ids, dijkstra.unreached_ids());
    // correctness of path from source to each target
    for (int64_t target_id = 0; target_id < n_vertices; ++target_id)
    {
        if (target_id == source_id)
            continue;
        cout << "finding all paths from source " << source_id << " to target " << target_id << endl;
        if (target_all_paths[target_id].empty())
            cout << "no path found" << endl;
        else
        {
            EXPECT_TRUE(reached_ids.find(target_id) != reached_ids.end());
            auto shortest_path_info = ShortestOfAllPaths<typename DijkstraType::GraphType>(test_graph, target_all_paths[target_id]);
            EXPECT_EQ(shortest_path_info.first, test_graph.Dsrc(target_id));
            cout << "shortest distance = " << test_graph.Dsrc(target_id) << ", ";
            int64_t current_id = target_id;
            for (auto rit = shortest_path_info.second.crbegin();
                 rit != shortest_path_info.second.crend(); ++rit)
            {
                EXPECT_EQ(*rit, current_id);
                current_id = test_graph.ParentVertexId(current_id);
            }
            cout << "shortest path = " << mcp3d::JoinVector(shortest_path_info.second) << endl;
        }
    }
}

TEST(MDijkstra, FindPathMultiTargets)
{
    cout << "executing MDijkstraFindPathMultiTargetsImpl 5 times" << endl;
    for (int i = 0; i < 5; ++i)
        MDijkstraFindPathMultiTargetsImpl();
}

void MDijkstraUpdateBoundaryImpl()
{
    using DijkstraType = mcp3d::MDijkstra<TestDoubleGraph>;
    DijkstraType dijkstra;
    typename DijkstraType::GraphType test_graph;
    int n_vertices = 20;
    double n_edge_percent = 0.15;
    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    auto n_edge = (int)(n_edge_percent * n_vertices * n_vertices);
    cout << "random test case: " << n_vertices << " vertices, " << n_edge << " edges" << endl;
    test_graph.GenerateRandomGraph(n_vertices, n_edge);
    // take random source id
    uniform_int_distribution<int64_t> vertex_distribution((int64_t)0, (int64_t)(n_vertices - 1));
    int64_t source_id = vertex_distribution(generator);
    // find shortest paths between all target and source
    dijkstra.FindPaths(test_graph, source_id, {}, true);

    unordered_map<int64_t, bool> has_children;
    // for each target_id, find all paths between source and target
    unordered_map<int64_t, vector<VertexPath>> target_all_paths;
    for (int64_t target_id = 0; target_id < n_vertices; ++target_id)
    {
        // vertices unreachable from source are not considered as leaf vertices
        if (dijkstra.unreached_ids().find(target_id) != dijkstra.unreached_ids().end())
            continue;
        // initialize has_children entry as false
        if (has_children.find(target_id) == has_children.end())
            has_children[target_id] =false;
        // for target id not equal to source id, find all path between source and target
        if (target_id == source_id)
            continue;
        target_all_paths[target_id] = test_graph.AllPaths(source_id, target_id);
        // any vertex that appeared on a path to the target vertex has at least
        // one child
        auto shortest_path_info = ShortestOfAllPaths<typename DijkstraType::GraphType>(test_graph, target_all_paths[target_id]);
        for (size_t i = 0; i < shortest_path_info.second.size() - 1; ++i)
            has_children[shortest_path_info.second[i]] = true;
    }
    size_t n_leaves = 0;
    for (const auto& vertex_has_children: has_children)
    {
        if (vertex_has_children.second)
            continue;
        ++n_leaves;
        EXPECT_TRUE(dijkstra.boundary().find(vertex_has_children.first) != dijkstra.boundary().end());
    }
    EXPECT_EQ(n_leaves, dijkstra.boundary().size());
    cout << n_leaves << " vertices in dijkstra boundary" << endl;
}

TEST(MDijkstra, UpdateBoundary)
{
    cout << "executing MDijkstraUpdateBoundaryImpl 5 times" << endl;
    for (int i = 0; i < 5; ++i)
        MDijkstraUpdateBoundaryImpl();
}


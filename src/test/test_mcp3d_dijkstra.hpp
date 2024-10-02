//
// Created by muyezhu on 10/6/19.
//

#ifndef MCP3D_TEST_MCP3D_DIJKSTRA_HPP
#define MCP3D_TEST_MCP3D_DIJKSTRA_HPP

#include <fstream>
#include <string>
#include <limits>
#include <chrono>
#include <random>
#include <vector>
#include <unordered_map>
#include "common/mcp3d_macros.hpp"
#include "algorithm/mcp3d_algorithm_macros.hpp"


namespace mcp3d
{

namespace test
{

/// vertex class used for unit test
template <typename EdgeCostType>
class TestDijkstraVertex: public mcp3d::MVertexBase<EdgeCostType>
{
public:
    TestDijkstraVertex(): edges_(std::unordered_map<int64_t, EdgeCostType>{}),
                          id_(-1), parent_id_(std::numeric_limits<int64_t>::lowest()),
                          dsrc_(std::numeric_limits<EdgeCostType>::infinity()), is_frozen_(false) {}

    explicit TestDijkstraVertex(int64_t id): edges_(std::unordered_map<int64_t, EdgeCostType>{}),
                                             id_(id), parent_id_(std::numeric_limits<int64_t>::lowest()),
                                             dsrc_(std::numeric_limits<EdgeCostType>::infinity()), is_frozen_(false) {}

    bool StoreEdges() const override { return true; }

    bool operator==(const TestDijkstraVertex& other) const
    { return  id_ == other.id_ && parent_id_ == other.parent_id_ && edges_ == other.edges_; }

    bool operator<(const TestDijkstraVertex& other) const
    { return dsrc_ < other.dsrc_; }

    bool HasNeighbor(int64_t neighbor_id) const  { return edges_.find(neighbor_id) != edges_.end(); }

    int NOutEdges() const  { return (int)edges_.size(); }

    void RemoveNeighbor(int64_t neighbor_id)  { edges_.erase(neighbor_id); }

    // neighbor_id can be pre-existing or new
    void UpdateNeighbor(int64_t neighbor_id, EdgeCostType edge_cost)  { edges_[neighbor_id] = edge_cost; }

    int64_t id() const override { return id_; }

    int64_t parent_id() const  { return parent_id_; }

    EdgeCostType dsrc() const  { return dsrc_; }

    bool is_frozen() const  { return is_frozen_; }

    std::unordered_map<int64_t, EdgeCostType> edges() override  { return edges_; };

    std::unordered_map<int64_t, EdgeCostType>& edges_ref()   { return edges_; };

    const std::unordered_map<int64_t, EdgeCostType>& edges_ref() const  { return edges_; };

    void set_id(int64_t id) override  { id_ = id; }

    void set_parent_id(int64_t parent_id)  { parent_id_ = parent_id; }

    void set_dsrc(EdgeCostType dsrc)  { dsrc_ = dsrc; }

    void set_is_frozen(bool is_frozen)  { is_frozen_ = is_frozen; }
private:
    std::unordered_map<int64_t, EdgeCostType> edges_;
    int64_t id_, parent_id_;
    EdgeCostType dsrc_;
    bool is_frozen_;
};


class VertexPath
{
public:
    VertexPath(): vertex_sequence_(std::vector<int64_t>{}), vertex_ids_(std::unordered_set<int64_t>{}) {}

    VertexPath(const VertexPath& other): vertex_sequence_(other.vertex_sequence_), vertex_ids_(other.vertex_ids_) {}

    bool Empty() const  { return vertex_sequence_.empty(); }

    // number of unique vertices
    size_t NVertices() const  { return vertex_ids_.size(); }

    size_t NEdges() const  { return std::max((size_t)0, vertex_sequence_.size() - 1); }

    bool AppendVertexId(int64_t vertex_id, bool allow_ancestor = false)
    {
        EXPECT_TRUE(vertex_id >= 0);
        if (!allow_ancestor &&
            vertex_ids_.find(vertex_id) != vertex_ids_.end())
            return false;
        vertex_sequence_.push_back(vertex_id);
        vertex_ids_.insert(vertex_id);
        return true;
    }

    int64_t PathEndVertexId() const  { return Empty() ? -1 : vertex_sequence_.back(); }

    std::vector<int64_t> vertex_sequence() const  { return vertex_sequence_; }

    void PrintPath() const
    {
        std::cout << "printing VertexPath: ";
        for (int64_t vertex_id: vertex_sequence_)
            std::cout << vertex_id << " ";
        std::cout << std::endl;
    }

private:
    std::vector<int64_t> vertex_sequence_;
    std::unordered_set<int64_t> vertex_ids_;
};

/// graph class used for unit test
template <typename Vertex>
class TestDijkstraGraph: public mcp3d::MGraph<Vertex>
{
public:
    MCP3D_MGRAPH_ALIASES

    TestDijkstraGraph() = default;

    explicit TestDijkstraGraph(const std::string& vertex_file_path)
    {
        EXPECT_TRUE(mcp3d::IsFile(vertex_file_path));
        std::string line;
        std::ifstream in_file(vertex_file_path.c_str(), std::ifstream::in);
        std::vector<std::string> tokens;
        while (getline(in_file, line))
        {
            tokens = mcp3d::SplitString(mcp3d::Strip(line), ",");
            if (tokens.empty())
                continue;
            EXPECT_TRUE(tokens.size() % 2 == 1);
            size_t n_neighbors = (tokens.size() - 1) / 2;
            int64_t vertex_id = stol(tokens[0]);
            EXPECT_TRUE(Vertices().find(vertex_id) == Vertices().end());
            operator[](vertex_id) = TestDijkstraVertex<EdgeType>(vertex_id);
            for (size_t i = 0; i < n_neighbors; ++i)
            {
                int64_t neighbor_id = stol(tokens[i * 2 + 1]);
                double edge_cost = stod(tokens[i * 2 + 2]);
                EXPECT_FALSE(At(vertex_id).HasNeighbor(neighbor_id));
                At(vertex_id).UpdateNeighbor(neighbor_id, edge_cost);
                EXPECT_EQ(edge_cost, At(vertex_id).edges_ref()[neighbor_id]);
            }
            EXPECT_EQ((int)n_neighbors, At(vertex_id).NOutEdges());
        }
    }

    void Clear() override { Vertices().clear(); }

    // generate random graph with n_vertices vertices and n_edges edges
    // if overwrite_existing_edge is false, each edge updated must be new,
    // new candidate edge is regenerated until the candidate does not previously
    // exist in the graph
    void GenerateRandomGraph(int n_vertice, int n_edge, bool overwrite_existing_edge = false)
    {
        MCP3D_ASSERT(n_vertice > 0 && n_edge >= 0);
        Clear();
        // create n_vertice vertices
        for (int i = 0; i < n_vertice; ++i)
        {
            Vertices().emplace((int64_t)i, Vertex{});
            At((int64_t)i).set_id((int64_t)i);
        }
        int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        using EdgeType = typename Vertex::EdgeType;
        std::uniform_int_distribution<int64_t> vertex_distribution((int64_t)0, (int64_t)(n_vertice - 1));
        std::uniform_real_distribution<EdgeType> edge_distribution((EdgeType)1, (EdgeType)100);
        int64_t start_id, end_id;
        for (int i = 0; i < n_edge; ++i)
        {
            do
            {
                start_id = vertex_distribution(generator);
                end_id = vertex_distribution(generator);
                if (overwrite_existing_edge)
                    break;
            } while (At(start_id).edges_ref().find(end_id) != At(start_id).edges_ref().end());
            EdgeType edge_cost = edge_distribution(generator);
            At(start_id).UpdateNeighbor(end_id, edge_cost);
        }
    }

    bool IsDirected() const override { return true; }

    std::vector<int64_t> NeighborVertexIds(int64_t vertex_id) const override
    {
        std::vector<int64_t> neighbor_ids;
        for (const auto& edge: At(vertex_id).edges_ref())
            neighbor_ids.push_back(edge.first);
        return neighbor_ids;
    }

    void NonFrozenNeighborVertexIds(int64_t vertex_id, std::unique_ptr<int64_t []>& non_frozen_neighbor_ids) const
    {
        non_frozen_neighbor_ids = std::make_unique<int64_t []>(27);
        int n_non_frozen_neighbors = 0;
        for (const auto& edge: At(vertex_id).edges_ref())
            if (!IsFrozen(edge.first))
                non_frozen_neighbor_ids[++n_non_frozen_neighbors] = edge.first;
        non_frozen_neighbor_ids[0] = n_non_frozen_neighbors;
    }

    EdgeType EdgeCost(int64_t vertex_id0, int64_t vertex_id1) override { return At(vertex_id0).edges_ref().at(vertex_id1); }

    EdgeType EdgeCost(int64_t vertex_id0, int64_t vertex_id1) const override { return At(vertex_id0).edges_ref().at(vertex_id1); }

    std::unordered_map<int64_t, EdgeType> Edges(int64_t vertex_id) override { return At(vertex_id).edges(); };

    int64_t ParentVertexId(int64_t vertex_id) const  { return At(vertex_id).parent_id(); }

    EdgeType Dsrc(int64_t vertex_id) const { return At(vertex_id).dsrc(); }

    bool IsFrozen(int64_t vertex_id) const  { return At(vertex_id).is_frozen(); }

    void SetDsrc(int64_t vertex_id, EdgeType dsrc)  { At(vertex_id).set_dsrc(dsrc); }

    void SetParentVertexId(int64_t vertex_id, int64_t parent_vertex_id)  { At(vertex_id).set_parent_id(parent_vertex_id); }

    void SetIsFrozen(int64_t vertex_id, bool is_frozen)  { At(vertex_id).set_is_frozen(is_frozen); }

    size_t NEdges() const
    {
        size_t n_edges = 0;
        for (const auto& vertex: Vertices())
            n_edges += vertex.second.NOutEdges();
        return n_edges;
    }

    // find all path between src_id and target_id
    // the path will be src_id, v0, v1, ..., vn, target_id, this is true
    // when src_id = target_id as well
    std::vector<VertexPath> AllPaths(int64_t src_id, int64_t target_id) const
    {
        EXPECT_TRUE(HasVertex(src_id) && HasVertex(target_id));
        INIT_TIMER(0)
        TIC(0)
        std::vector<VertexPath> complete_paths, growing_paths;
        // initialize the src_id as start point for all paths
        VertexPath initial_path {};
        initial_path.AppendVertexId(src_id);
        growing_paths.push_back(initial_path);
        int64_t n_examined_paths = 0;
        while (!growing_paths.empty())
        {
            int64_t path_end_id = growing_paths[0].PathEndVertexId();
            // if a growing path has reached target_id, create complete paths
            // when source id = target id, the shortest path should be
            // [source id, target id]
            if (path_end_id == target_id && growing_paths[0].NEdges() >= 1)
                complete_paths.push_back(growing_paths[0]);
                // grow path by 1 edge if no cycle is detected. add the newly
                // grown path to the vector
            else
            {
                std::vector<int64_t> neighbor_ids(NeighborVertexIds(path_end_id));
                for (int64_t neighbor_id: neighbor_ids)
                {
                    VertexPath extend_path(growing_paths[0]);
                    // in the case when target id = source id,
                    // the target id is allowed to be appended into the path
                    // though its been seen before (as source id)
                    bool no_cycle = extend_path.AppendVertexId(neighbor_id, neighbor_id == target_id);
                    if (no_cycle)
                        growing_paths.push_back(extend_path);
                }
            }
            // erase the path at the very beginning of the vector.
            growing_paths.erase(growing_paths.begin());
            ++n_examined_paths;
        }
        TOC(0)
        std::cout << "number of growing paths enumerated = " << n_examined_paths << std::endl;
        REPORT_TIME_TO_COMPLETION("time to enumerate growing paths", 0)
        return complete_paths;
    }
};

}

}


#endif //MCP3D_TEST_MCP3D_DIJKSTRA_HPP

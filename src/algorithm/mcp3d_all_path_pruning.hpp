//
// Created by muyezhu on 9/19/19.
//

#ifndef MCP3D_MCP3D_ALL_PATH_PRUNING_GENERAL_HPP
#define MCP3D_MCP3D_ALL_PATH_PRUNING_GENERAL_HPP


#include <cstdint>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include "common/mcp3d_macros.hpp"
#include "image/mcp3d_image_utils.hpp"
#include "image/mcp3d_image_maths.hpp"
#include "mcp3d_vertex.hpp"
#include "mcp3d_graph.hpp"
#include "mcp3d_dijkstra.hpp"


namespace mcp3d
{

template <typename Graph>
class MHierarchicalSegment
{
public:
    using EdgeType = typename Graph::EdgeType;

    MHierarchicalSegment(): parent_segment_(nullptr),
                           leaf_id_(std::numeric_limits<int64_t>::lowest()),
                           root_id_(std::numeric_limits<int64_t>::lowest()),
                           length_(0.0), level_(1) {}

    MHierarchicalSegment(int64_t leaf_id, int64_t root_id, double length, int level):
            parent_segment_(nullptr), leaf_id_(leaf_id), root_id_(root_id),
            length_(length), level_(level) {}

    void VertexSequence(Graph& graph, std::vector<int64_t> & vertex_sequence)
    {
        if(leaf_id_ < 0 || root_id_ < 0)
            return;
        int64_t vertex_id = leaf_id_;
        while(vertex_id != root_id_)
        {
            vertex_sequence.push_back(vertex_id);
            vertex_id = graph.ParentVertex(vertex_id);
        }
        vertex_sequence.push_back(root_id_);
    }

    int64_t leaf_id() const  { return leaf_id_; }

    int64_t root_id() const  { return root_id_; }

    EdgeType length() const  { return length_; }

    void set_parent_segment(MHierarchicalSegment* parent_segment)  { parent_segment_ = parent_segment; }

private:
    MHierarchicalSegment* parent_segment_;
    int64_t leaf_id_, root_id_;   // its parent marker is in current segment's parent segment
    EdgeType length_;             // the length from leaf to root
    int level_;                   // the segments number from leaf to root
};


/// performs routine fast marching
/// additionally maintains the narrow band fronts. each vertex retrieved
/// from the heap is inserted into the front, while its parent is removed
/// propose: at dijkstra end, remove frontier nodes that are not border voxels
/// from frozen, and run one more step of dijkstra from each non border frontier
/// voxel
template <typename Vertex, typename VType = uint16_t>
class MAllPathPruning
{
public:
    using VertexType = Vertex;
    using EdgeType = typename Vertex::EdgeType;
    using GraphType = MGraphAllPath<Vertex, VType>;
    using HierarchicalSegmentType = MHierarchicalSegment<GraphType>;

    MAllPathPruning(): graph_(GraphType{}), dijkstra_(MDijkstra<GraphType>{})  {}

    void Clear();

    // placeholder for input transformation functions. should call
    // generic transformation functions defined else where
    void Transform(void* input, const std::string& method = "distance");

    // actual "travel time" impelemnetation is not solution to the quadratic
    // formulation, but based on graph edge costs. therefroe this is in essence
    // dijkstra's algothrim
    // background voxels in input are expected to be 0
    bool ConstructPaths(void* input, int input_zdim, int input_ydim, int input_xdim, int srcz, int srcy, int srcx, int cnn_type = 3);

    bool ConstructPaths(void* input, int input_zdim, int input_ydim, int input_xdim, double srcz, double srcy, double srcx, int cnn_type = 3);

    void HierarchicalPrune(EdgeType length_thresh = -1.0, const std::string& distance_metric = "intensity");

    GraphType& Graph()  { return graph_; }

private:
    void ConstructHierarchicalSegments(const std::string& distance_metric);

    void PruneShortSegments(EdgeType length_thresh);

    void ComputeVertexRadii();

    GraphType graph_;
    MDijkstra<GraphType> dijkstra_;
    /// leaf_id: hierarchical_segment with leaf_id as its leaf_id_ member field
    std::unordered_map<int64_t, HierarchicalSegmentType> hierarchical_segments_;
};

}

template <typename Vertex, typename VType>
void mcp3d::MAllPathPruning<Vertex, VType>::Clear()
{
    graph_.Clear();
    dijkstra_.Clear();
};

template <typename Vertex, typename VType>
bool mcp3d::MAllPathPruning<Vertex, VType>::ConstructPaths(void* input, int input_zdim, int input_ydim, int input_xdim,
                                                           double srcz, double srcy, double srcx, int cnn_type)
{

    return ConstructPaths(input, input_zdim, input_ydim, input_xdim,
                          (int)std::round(srcz), (int)std::round(srcy),
                          (int)std::round(srcx), cnn_type);
}

template <typename Vertex, typename VType>
bool mcp3d::MAllPathPruning<Vertex, VType>::ConstructPaths(void* input, int input_zdim, int input_ydim, int input_xdim,
                                                           int srcz, int srcy, int srcx, int cnn_type)
{
    INIT_TIMER(0)
    TIC(0)
    graph_.Init(input, input_zdim, input_ydim, input_xdim, srcz, srcy, srcx, cnn_type);
    TOC(0)
    REPORT_TIME_TO_COMPLETION("MAllPathPruning graph initialization", 0)
    // dijkstra with boundary tracking
    bool success = dijkstra_.FindPaths(graph_, graph_.SrcId(), std::unordered_set<int64_t>{}, false);
    if (!success)
        MCP3D_MESSAGE("not all foreground voxels are reachable from source node")
}

template <typename Vertex, typename VType>
void mcp3d::MAllPathPruning<Vertex, VType>::HierarchicalPrune(EdgeType length_thresh, const std::string& distance_metric)
{
    ConstructHierarchicalSegments(distance_metric);

    size_t n_segments_preprune = hierarchical_segments_.size();
    PruneShortSegments(length_thresh);
    size_t n_segments_postprune = hierarchical_segments_.size();
    std::cout << "pruned by length_thresh (segment number) : "
              << n_segments_preprune << " - " << n_segments_preprune - n_segments_postprune
              << " = " << n_segments_postprune << std::endl;

    std::multimap<EdgeType, HierarchicalSegmentType*> seg_length_map;
    for(auto& item: hierarchical_segments_)
        seg_length_map.emplace(item.second.length(), &(item.second));


}

template <typename Vertex, typename VType>
void mcp3d::MAllPathPruning<Vertex, VType>::ConstructHierarchicalSegments(const std::string& distance_metric)
{
    MCP3D_ASSERT(dijkstra_.find_path_complete())
    std::cout << "Construct hierarchical segments" << std::endl;
    int64_t src_id = graph_.SrcId();
    std::unordered_set<int64_t> reached_ids(dijkstra_.ReachedTargetIds());
    std::unordered_map<int64_t, int64_t>& dijkstra_boundary = dijkstra_.boundary();
    std::unordered_map<int64_t, int64_t> n_children,       // number of children for each vertex derived from dijkstra's
                                         vertex_leaf_ids;  // the leaf id of each vertex, determined by argmax_leaf(path_weights(vertex, leaf))
    std::unordered_map<int64_t, EdgeType> path_weights;     // weight of the path vertex -> leaf

    // initialize variables
    n_children[src_id] = 0;
    vertex_leaf_ids[src_id] = std::numeric_limits<int64_t>::lowest();
    path_weights[src_id] = std::numeric_limits<VType>::lowest();
    for (const auto& reached_id: reached_ids)
    {
        n_children[reached_id] = 0;
        vertex_leaf_ids[reached_id] = std::numeric_limits<int64_t>::lowest();
        path_weights[reached_id] = std::numeric_limits<VType>::lowest();
    }

    // for each vertex, obtain number of child vertices
    for (const auto& reached_id: reached_ids)
    {
        int64_t parent_id = graph_.ParentVertexId(reached_id);
        ++n_children.at(parent_id);
    }
    // remove parent (std::numeric_limits<int64_t>::lowest()) of source vertex
    n_children.erase(std::numeric_limits<int64_t>::lowest());
    vertex_leaf_ids.erase(std::numeric_limits<int64_t>::lowest());
    path_weights.erase(std::numeric_limits<int64_t>::lowest());

    // assign a leaf vertex to each vertex, and compute the weights of the path
    // from leaf to the vertex using intensity values
    EdgeType candidate_weight;
    // vertices in dijkstra boundary set are also the leaf vertices
    for (const auto& item: dijkstra_boundary)
    {
        int64_t leaf_id = item.first, vertex_id, parent_id;
        vertex_leaf_ids.at(leaf_id) = leaf_id;
        path_weights.at(leaf_id) = graph_.NormalizedDataValue(leaf_id);
        vertex_id = leaf_id;
        while (graph_.HasParent(vertex_id))
        {
            parent_id = graph_.ParentVertexId(vertex_id);
            candidate_weight = path_weights.at(vertex_id) + graph_.NormalizedDataValue(parent_id);
            if (candidate_weight >= path_weights.at(parent_id))
            {
                path_weights.at(parent_id) = candidate_weight;
                vertex_leaf_ids.at(parent_id) = leaf_id;
            }
            vertex_id = parent_id;
        }
    }

    // create hierarchical segments
    for (const auto& item: dijkstra_boundary)
    {
        int64_t leaf_id = item.first;
        int level = 1;
        // start_id: start id in segment. the segment ends with leaf id
        int64_t root_id = leaf_id, start_parent_id;
        while (graph_.HasParent(root_id))
        {
            start_parent_id = graph_.ParentVertexId(root_id);
            if (vertex_leaf_ids.at(start_parent_id) != leaf_id)
                break;
            if (n_children.at(root_id) >= 2)
                ++level;
            root_id = start_parent_id;
        }
        hierarchical_segments_.emplace(leaf_id, HierarchicalSegmentType(leaf_id, root_id, path_weights[root_id], level));
    }
    // enumerate parent segment of each hierarchical segment
    for (auto& item: hierarchical_segments_)
    {
        HierarchicalSegmentType& segment = item.second;
        int64_t root_id = segment.root_id();
        if (graph_.HasVertex(root_id))
        {
            int64_t root_parent_id = graph_.ParentVertexId(root_id);
            int64_t root_parent_leaf_id = vertex_leaf_ids.at(root_parent_id);
            segment.set_parent_segment(&(hierarchical_segments_.at(root_parent_leaf_id)));
        }
    }
}

template <typename Vertex, typename VType>
void mcp3d::MAllPathPruning<Vertex, VType>::PruneShortSegments(EdgeType length_thresh)
{
    if (length_thresh < 0)
        return;
    std::vector<int64_t> short_segment_ids;
    for (const auto& item: hierarchical_segments_)
        if (item.second.length() < length_thresh)
            short_segment_ids.push_back(item.first);
    for (const auto& short_segment_id: short_segment_ids)
        hierarchical_segments_.erase(short_segment_id);
}

template <typename Vertex, typename VType>
void mcp3d::MAllPathPruning<Vertex, VType>::ComputeVertexRadii()
{
    VType data_ptr = graph_.data();
    for (auto& segments_: hierarchical_segments_)
    {

    }
}


#endif //MCP3D_MCP3D_ALL_PATH_PRUNING_GENERAL_HPP

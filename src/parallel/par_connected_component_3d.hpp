//
// Created by muyezhu on 12/11/17.
//

#ifndef MCP3D_PAR_CONNECTED_COMPONENT_3D_HPP
#define MCP3D_PAR_CONNECTED_COMPONENT_3D_HPP

#include <vector>
#include <unordered_map>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "algorithm/mcp3d_connected_component_3d.hpp"

static const int MPI_CC3D_W = 8192;
static const int MPI_CC3D_H = 128;
static const int MPI_CC3D_Z = 100;
namespace mcp3d
{

class ComponentRecord
{
public:
    explicit ComponentRecord(int32_t l);
    void AddPoint(Position3D<int32_t>& p) { AddPoint(std::move(p)); }
    void AddPoint(Position3D<int32_t>&& p);
    void CombineWith(ComponentRecord &other);
    void Save(std::string out_dir);
    void ComputeCenter();
    bool operator== (const ComponentRecord &other);
    int32_t label() { return label_; };
    int32_t mass() { return mass_; };
    Position3D<double> center_of_mass() { return center_of_mass_; };
    std::vector<Position3D<int32_t>>& coordinates() { return coordinates_; }
private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &ar, const uint32_t version)
    {
        ar & label_;
        ar & mass_;
        ar & center_of_mass_;
        ar & coordinates_;
    }
    int32_t label_, mass_;
    Position3D<double> center_of_mass_;
    std::vector<Position3D<int32_t>> coordinates_;
};

/// distributed connected component label. a large image I of dimensions
/// (Z, Y, X) is divided into blocks of dimensions (Z, img_h, img_w), with
/// img_h and img_w being user supplied values. to be compliant with tiff image
/// tile size requirement, both img_h and img_w should be multiples of 16.
/// the current code only implements neighbor face communication, results will
/// only be correct when no edge or corner communication is required. therefore
/// its important to use blocks with dimensions (Z, img_h, X), so each process
/// only has face neighbor in y axis direction(or Z, Y, img_w is correct too)
class ParCCLabel3D
{
public:
    explicit ParCCLabel3D(const std::string& tissue_dir,
                          int img_w = MPI_CC3D_W, int img_h = MPI_CC3D_H,
                          int foreground = 1);
    /// the z level images are read as first axis, so [0, :, :] is first image
    /// the axis are in z, y, x order. [i, j, k] address voxel at Zi, Yj, Xk
    /// y: image height direction, x: image width direction
    /// the block has a one pixel wide boundary area whenever available. the
    /// boundary pixels are used to resolve label equivalence across processors
    /// the block with block coordinate (bi, bj, bk) will read the image volume
    /// I[:, bi * img_h - 1: (bi + 1) * img_h, bj * img_w - 1, (bj + 1) * img_w]
    void ReadInput();
    void SetForeground();
    void LabelLocal();
    void GlobalizeLabel();
    /// each processor communicates with its neighbors the labels of pixels with
    /// processor boundary region. only face neighbor communication implemented.
    /// both 1 and -1 directions needed to be executed for correct processor
    /// boundary label resolution. the label equivalence map is used in next
    /// stage communication
    void CommunicateNeighborLabels();
    /// merge label equivalence map from each processor in reduction tree
    /// pattern. the root broadcasts full map
    void GlobalLabelReduction();
    void FinalizeLabel();
    int32_t CountGlobalCC();
    void WriteLabeledImages();
    void WriteComponents(std::string out_dir);
    void QuantifyComponenets();
    int process_id()  { return process_id_; }
    int n_proc() { return n_proc_; }
    std::vector<int64_t> tissue_lens()
    { return std::vector<int64_t>({tissue_zlen_, tissue_ylen_, tissue_xlen_}); }
    std::vector<int64_t>& block_dims()  { return block_dims_; }
    ~ParCCLabel3D() { delete[] label; }
private:
    int ProcessIDFromXYZID(int64_t target_xid, int64_t target_yid, int64_t target_zid);
    void BlockNonPadSurfacesIndices(int32_t *i_north, int32_t *i_south,
                                    int32_t *j_west, int32_t *j_east,
                                    int32_t *k_front, int32_t *k_back);
    void BlockSurfaceIndices(int32_t* i_north, int32_t* i_south,
                             int32_t* j_west, int32_t* j_east,
                             int32_t* k_front, int32_t* k_back);
    /// direction: 1 or -1, message sent to positive vs negative direction of axis
    /// e.g. axis = 0, direction = -1: (z, y, x) sending to (z - 1, y, x)
    /// the top left pixel of first image in the z stack is coordinate origin
    void CommunicateNeighborFaceLabels(int axis, int direction);
    /// send_direction is direction of the sent message from the sender
    void ResolveNeighborFaceLabels(int32_t* rec_buf, int axis,
                                   int send_direction);
    void MergeGlobalLabels(int n_label_pairs, int32_t* rec_buf);
    int64_t BlockVecCoorToScalar(int32_t i, int32_t j, int32_t k);
    /// only traverse non padding voxels
    void BuildLocalComponentRecords();
    void WriteLocalComponents(std::string out_dir);
    void CommunicateComponentRecords();
    bool is_worker_, label_finalized_, is_master_recorder_;
    std::string tissue_dir_;
    std::vector<std::string> img_paths_;
    int process_id_, n_proc_, required_worker_num_, required_total_num_;
    int64_t block_xid_, block_yid_, block_zid_;
    int64_t tissue_xlen_, tissue_ylen_, tissue_zlen_;
    int64_t block_xlen_, block_ylen_, block_zlen_;
    int32_t tile_len_, total_tile_num_;
    int64_t tissue_xdim_, tissue_ydim_, tissue_zdim_;
    int64_t block_xpad_, block_ypad_, block_zpad_;
    int64_t block_total_voxels_;
    /// (block_zlen_ + block_zpad_, block_ylen_ + block_ypad_,
    /// block_xlen_ + block_xpad_)
    std::vector<int64_t> block_dims_;
    int32_t local_n_cc_, global_n_cc_, global_label_offset_;
    int foreground_;
    int32_t* label;
    uint8_t* input;
    std::unordered_map<int32_t, int32_t> global_eq_labels_;
    /// border_labels_: only useful before global label reduction. set of local
    /// labels falling in border regions
    /// final_labels_: only useful after global label reduction. set of labels
    /// in each processor's block after label reduction
    std::unordered_set<int32_t> border_labels_, final_labels_;
    std::unordered_map<int32_t, ComponentRecord> component_records_;
};
}


#endif //MCP3D_PAR_CONNECTED_COMPONENT_3D_HPP

//
// Created by muyezhu on 12/11/17.
//
#include <iostream>
#include <fstream>
#include <chrono>
#include <new>
#include <utility>
#include <cstring>
#include <cmath>
#include <unordered_map>
#include <mpi.h>
#include <boost/serialization/vector.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include "image_interface/mcp3d_tiff_utils.hpp"
#include "par_connected_component_3d.hpp"

using namespace std;

enum NeighborBlockType
{
    FACE, EDGE, CORNER
};

mcp3d::ComponentRecord::ComponentRecord(int32_t l):
        label_(l), mass_(1), center_of_mass_(Position3D<double>())
{
    if (label_ < 0 ) MCP3D_INVALID_ARGUMENT("component label must be positive")
}

bool mcp3d::ComponentRecord::operator==(const ComponentRecord &other)
{
    if (label_ != other.label_)
    {
        cout << "components have different labels" << endl;
        return false;
    }
    if (center_of_mass_ != other.center_of_mass_)
    {
        cout << "components have different center of masses" << endl;
        return false;
    }
}

void mcp3d::ComponentRecord::AddPoint(Position3D<int32_t>&& p)
{
    // move constructor
    coordinates_.push_back(p);
    ++mass_;
}

void mcp3d::ComponentRecord::CombineWith(ComponentRecord &other)
{
    if (label_ != other.label()) MCP3D_INVALID_ARGUMENT(
            "can not combine components with non identical labels")
    mass_ += other.mass();
    coordinates_.insert(coordinates_.end(), other.coordinates().begin(), other.coordinates().end());
}

void mcp3d::ComponentRecord::Save(string out_dir)
{
    string record_path = out_dir + "/component_" + to_string(label_);
    ofstream ofs(record_path);
    boost::archive::binary_oarchive boa(ofs);
    boa << *this;
}

// TODO: use MImage instance
mcp3d::ParCCLabel3D::ParCCLabel3D(const std::string &tissue_dir,
                                  int img_w, int img_h, int foreground):
              is_worker_(false), label_finalized_(false), is_master_recorder_(false),
              tissue_dir_(tissue_dir), block_xlen_(img_w), block_ylen_(img_h),
              local_n_cc_(0), global_label_offset_(0), foreground_(foreground),
              label(nullptr), input(nullptr)
{
    if (!mcp3d::IsDir(tissue_dir)) MCP3D_INVALID_ARGUMENT_FINALIZE(
            "tissue dir doesn't exist: " + tissue_dir)
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized == 0) MCP3D_RUNTIME_ERROR_FINALIZE(
            "ParCCLabel3D can not be instantiated without prior call to MPI_Init")
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc_);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id_);
    if (block_xlen_ % 16 != 0 && block_ylen_ % 16 != 0) MCP3D_INVALID_ARGUMENT_FINALIZE(
            "at least one of block_xlen and block_ylen must be multiples of 16")
    tile_len_ = (int32_t)min((int64_t)mcp3d::TIFFTILE_XDIM,
                             min(block_xlen_, block_ylen_));
    total_tile_num_ = (tissue_xlen_ / tile_len_ +
                       tissue_xlen_ % tile_len_ != 0 ? 1 : 0) *
                      (tissue_ylen_ / tile_len_ +
                       tissue_ylen_ % tile_len_ != 0 ? 1 : 0);
    mcp3d::TiffDirectoryInfo tiff_info(img_paths_[0]);
    tissue_xlen_ = (int)tiff_info.image_width,
    tissue_ylen_ = (int)tiff_info.image_height,
    tissue_zlen_ = (int)img_paths_.size();
    MCP3D_ASSERT(block_xlen_ <= tissue_xlen_ && block_ylen_ <= tissue_ylen_)
    block_zlen_ = (int)img_paths_.size();
    tissue_xdim_ = tissue_xlen_ % block_xlen_ == 0 ?
                  tissue_xlen_ / block_xlen_ : tissue_xlen_ / block_xlen_ + 1;
    tissue_ydim_ = tissue_ylen_ % block_ylen_ == 0 ?
                     tissue_ylen_ / block_ylen_ : tissue_ylen_ / block_ylen_ + 1;
    tissue_zdim_ = tissue_zlen_ % block_zlen_ == 0 ?
                  tissue_zlen_ / block_zlen_ : tissue_zlen_ / block_zlen_ + 1;
    required_worker_num_ = (int)(tissue_xdim_ * tissue_ydim_ * tissue_zdim_);
    // all workers plus one master recorder who holds no array data
    required_total_num_ = required_worker_num_ + 1;
    if (n_proc_ < required_total_num_)
        cout << "total proc num  = " << n_proc_
             << ", less than required total num = "
             << required_total_num_ << endl;
    if (process_id_ < required_worker_num_)
        is_worker_ = true;
    if (process_id_ == required_worker_num_)
        is_master_recorder_ = true;
    if (!is_worker_)
        return;
    block_zid_ = process_id_ / (tissue_xdim_ * tissue_ydim_);
    block_yid_ = process_id_ / tissue_xdim_ % tissue_ydim_;
    block_xid_ = process_id_ % tissue_xdim_;
    cout << "process id = " << process_id_ << ", zid = "
         << block_zid_ << ", yid = " << block_yid_
         << ", xid = " << block_xid_ << endl;
    if (tissue_xdim_ == 1)
        block_xpad_ = 0;
    else
    {
        if (block_xid_ == 0 || block_xid_ == tissue_xdim_ - 1)
            block_xpad_ = 1;
        else
            block_xpad_ = 2;
    }
    if (tissue_ydim_ == 1)
        block_ypad_ = 0;
    else
    {
        if (block_yid_ == 0 || block_yid_ == tissue_ydim_ - 1)
            block_ypad_ = 1;
        else
            block_ypad_ = 2;
    }
    if (tissue_zdim_ == 1)
        block_zpad_ = 0;
    else
    {
        if (block_zid_ == 0 || block_zid_ == tissue_zdim_ - 1)
            block_zpad_ = 1;
        else
            block_zpad_ = 2;
    }
    block_total_voxels_ = (block_xlen_ + block_xpad_) *
                          (block_ylen_ + block_ypad_) *
                          (block_zlen_ + block_zpad_);
    block_dims_.push_back(block_zlen_ + block_zpad_);
    block_dims_.push_back(block_ylen_ + block_ypad_);
    block_dims_.push_back(block_xlen_ + block_xpad_);
}

void mcp3d::ParCCLabel3D::ReadInput()
{
    if (!is_worker_)
        return;
    cv::Mat m;
    int64_t subimg_bytes = (block_xlen_ + block_xpad_) * (block_ylen_ + block_ypad_);
    input = new (nothrow) uint8_t[block_total_voxels_];
    if (!input) MCP3D_BAD_ALLOC("can not allocate memory for input array")
    for (size_t i = 0; i < img_paths_.size(); ++i)
    {
        mcp3d::SubImgRead(img_paths_[i], m,
                          block_xlen_ + block_xpad_, block_ylen_ + block_ypad_,
                          max(0l, block_xid_ * block_xlen_ - 1),
                          max(0l, block_yid_ * block_ylen_ - 1));
        memcpy(input + i * subimg_bytes, m.ptr(), size_t(subimg_bytes) * sizeof(uint8_t));
    }
}

void mcp3d::ParCCLabel3D::SetForeground()
{
    if (!is_worker_)
        return;
    if (foreground_ > 0)
        for (int64_t i = 0; i < block_total_voxels_; ++i)
            if (input[i] != foreground_)
                input[i] = 0;
}

void mcp3d::ParCCLabel3D::LabelLocal()
{
    if (!is_worker_)
        return;
    pair<int32_t, int32_t*> label_results = mcp3d::LabelCC3D<uint8_t, int32_t>(input, block_dims_);
    local_n_cc_ = label_results.first;
    label = label_results.second;
    cout << "process id = " << process_id_ << " local_ncc = " << local_n_cc_ << endl;
    delete[] input;
    MPI_Allreduce(&local_n_cc_, &global_n_cc_, 1, MPI_INT32_T, MPI_SUM, MPI_COMM_WORLD);
    cout << "process id = " << process_id_ << " global_ncc = " << global_n_cc_ << endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

void mcp3d::ParCCLabel3D::GlobalizeLabel()
{
    if (!is_worker_)
    {
        cout << "idle process, exiting from GlobalizeLabel()" << endl;
        return;
    }
    int tag = 1;
    MPI_Status status;
    int n_worker = min(required_worker_num_, n_proc_);
    if (process_id_ == 0)
        MPI_Send(&local_n_cc_, 1, MPI_INT32_T, process_id_ + 1, tag, MPI_COMM_WORLD);
    else if (process_id_ < n_worker - 1)
    {
        MPI_Recv(&global_label_offset_, 1, MPI_INT32_T, process_id_ - 1, tag, MPI_COMM_WORLD, &status);
        int32_t next_global_offset = global_label_offset_ + local_n_cc_;
        MPI_Send(&next_global_offset, 1, MPI_INT32_T, process_id_ + 1, tag, MPI_COMM_WORLD);
    }
    else if (process_id_ == n_worker - 1)
    {
        MPI_Recv(&global_label_offset_, 1, MPI_INT32_T, process_id_ - 1, tag, MPI_COMM_WORLD, &status);
    }
    if (is_worker_)
        for (int64_t i = 0; i < block_total_voxels_; ++i)
            if (label[i] > 0)
                label[i] += global_label_offset_;
    cout << "process id = " << process_id_ << ", global offset = " << global_label_offset_ << endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

void mcp3d::ParCCLabel3D::BlockSurfaceIndices(int32_t* i_north, int32_t* i_south,
                                              int32_t* j_west, int32_t* j_east,
                                              int32_t* k_front, int32_t* k_back)
{
    *i_north = 0;
    *j_west = 0;
    *k_front = 0;
    *i_south = (int32_t)block_zlen_ + (int32_t)block_zpad_ - 1;
    *j_east = (int32_t)block_ylen_ + (int32_t)block_ypad_ - 1;
    *k_back = (int32_t)block_xlen_ + (int32_t)block_xpad_ - 1;
}

void mcp3d::ParCCLabel3D::BlockNonPadSurfacesIndices(int32_t *i_north,
                                                     int32_t *i_south,
                                                     int32_t *j_west,
                                                     int32_t *j_east,
                                                     int32_t *k_front,
                                                     int32_t *k_back)
{
    if (tissue_zdim_ == 1)
    {
        *i_north = 0;
        *i_south = (int32_t)block_zlen_ - 1;
    }
    else
    {
        if (block_zid_ == 0)
        {
            *i_north = 0;
            *i_south = (int32_t)block_zlen_ - 1;
        }

        else
        {
            *i_north = 1;
            *i_south = (int32_t)block_zlen_;
        }
    }

    if (tissue_ydim_ == 1)
    {
        *j_west = 0;
        *j_east = (int32_t)block_ylen_ - 1;
    }
    else
    {
        if (block_yid_ == 0)
        {
            *j_west = 0;
            *j_east = (int32_t)block_ylen_ - 1;
        }

        else
        {
            *j_west = 1;
            *j_east = (int32_t)block_ylen_;
        }
    }
    if (tissue_xdim_ == 1)
    {
        *k_front = 0;
        *k_back = (int32_t)block_xlen_ - 1;
    }
    else
    {
        if (block_xid_ == 0)
        {
            *k_front = 0;
            *k_back = (int32_t)block_xlen_ - 1;
        }
        else
        {
            *k_front = 1;
            *k_back = (int32_t)block_xlen_ ;
        }
    }
}

int64_t mcp3d::ParCCLabel3D::BlockVecCoorToScalar(int32_t i, int32_t j, int32_t k)
{
    int64_t addr = (int64_t)i * block_dims_[1] * block_dims_[2] +
                 (int64_t)j * block_dims_[2] + (int64_t)k;
    return addr;
}

int mcp3d::ParCCLabel3D::ProcessIDFromXYZID(int64_t target_xid,
                                            int64_t target_yid,
                                            int64_t target_zid)
{
    if (target_xid < 0 || target_xid >= tissue_xdim_) MCP3D_INVALID_ARGUMENT(
            "target process xid out of range")
    if (target_yid < 0 || target_yid >= tissue_ydim_) MCP3D_INVALID_ARGUMENT(
            "target process yid out of range")
    if (target_zid < 0 || target_zid >= tissue_zdim_) MCP3D_INVALID_ARGUMENT(
            "target process zid out of range")
    return (int)(target_zid * tissue_xdim_ * tissue_ydim_ +
                 target_yid * tissue_xdim_ + target_xid);
}

void mcp3d::ParCCLabel3D::ResolveNeighborFaceLabels(int32_t* rec_buf, int axis,
                                                    int send_direction)
{
    if (!is_worker_)
        return;
    vector<pair<int32_t, int32_t>> nonpad_surface_ind({make_pair(-1, -1),
                                                   make_pair(-1, -1),
                                                   make_pair(-1, -1)});
    BlockNonPadSurfacesIndices(&nonpad_surface_ind[0].first,
                               &nonpad_surface_ind[0].second,
                               &nonpad_surface_ind[1].first,
                               &nonpad_surface_ind[1].second,
                               &nonpad_surface_ind[2].first,
                               &nonpad_surface_ind[2].second);
    int fixed_ind;
    if (send_direction == 1)
        fixed_ind = nonpad_surface_ind[axis].first;
    else
        fixed_ind = nonpad_surface_ind[axis].second;
    int n_cc = rec_buf[0];
    if (n_cc == CC_EMPTY)
    {
        cout << "process " << process_id_ << " receive empty neighbor labels. fix ind = " << fixed_ind << endl;
        return;
    }
    cout << "process " << process_id_ << " receive " << n_cc << " labels. fix ind = " << fixed_ind << endl;
    int32_t l_neighbor, l_self;
    int32_t i, j, k, c1, c2;
    int64_t addr;
    cout << "prrocess " << process_id_ << " l neighbor and l self: ";
    unordered_set<string> self_neighbor;
    for (int p = 1; p < 3 * n_cc + 1; p += 3)
    {
        l_neighbor = rec_buf[p];
        c1 = rec_buf[p + 1];
        c2 = rec_buf[p + 2];
        if (axis == 0)
        {
            i = fixed_ind;
            j = c1;
            k = c2;
        }
        else if (axis == 1)
        {
            i = c1;
            j = fixed_ind;
            k = c2;
        }
        else
        {
            i = c1;
            j = c2;
            k = fixed_ind;
        }
        addr = BlockVecCoorToScalar(i, j, k);
        l_self = label[addr];
        border_labels_.insert(l_self);
        if (border_labels_.find(l_neighbor) == border_labels_.end())
        {
            border_labels_.insert(l_neighbor);
            global_eq_labels_[l_neighbor] = l_self;
            global_eq_labels_[l_self] = CC_EMPTY;
        }
        else
        {   // parent of l_neighbor is not l_self, move l_self and its family to parent of l_neighbor
            if (global_eq_labels_[l_neighbor] != l_self)
            {
                int32_t par_lself = global_eq_labels_[l_self];
                if (par_lself == l_neighbor || par_lself == global_eq_labels_[l_neighbor])
                    continue;
                // if l_self is head, move l_self and its children to same root as l_neighbor
                if (par_lself == CC_EMPTY)
                {
                    for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
                    {
                        if (it->second == l_self)
                            it->second = global_eq_labels_[l_neighbor];
                    }
                    global_eq_labels_[l_self] = global_eq_labels_[l_neighbor];
                }
                // else move l_self and its parent and sibling to same root as l_neighbor
                else
                {
                    for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
                    {
                        if (it->second == par_lself)
                            it->second = global_eq_labels_[l_neighbor];
                    }
                    global_eq_labels_[par_lself] = global_eq_labels_[l_neighbor];
                    global_eq_labels_[l_self] = global_eq_labels_[l_neighbor];
                }
            }
        }
        if (self_neighbor.find(to_string(l_neighbor) + "_" + to_string(l_self)) == self_neighbor.end())
        {
            self_neighbor.insert(to_string(l_neighbor) + "_" + to_string(l_self));
            cout << "(" << l_neighbor << ", " << l_self <<")    ";
        }

    }
    cout << endl;
    cout << "process " << process_id_ << " receive from neighbor: ";
    for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
        cout << "(" << it->first << ", " << it->second << ")   ";
    cout << endl;
}

// local_l, j, k
void mcp3d::ParCCLabel3D::CommunicateNeighborFaceLabels(int axis, int direction)
{
    if (!is_worker_)
        return;
    if (axis != 0 && axis != 1 && axis != 2) MCP3D_INVALID_ARGUMENT(
            "axis out of range")
    if (direction != 1 && direction != -1) MCP3D_INVALID_ARGUMENT(
            "direction out of range")
    vector<int64_t> block_ids({block_zid_, block_yid_, block_xid_});
    vector<int64_t> block_lens({block_zlen_, block_ylen_, block_xlen_});
    vector<int64_t> tissue_dims({tissue_zdim_, tissue_ydim_, tissue_xdim_});
    if (tissue_dims[axis] == 1)
        return;
    vector<pair<int32_t, int32_t>> surface_ind({make_pair(-1, -1),
                                            make_pair(-1, -1),
                                            make_pair(-1, -1)});
    BlockSurfaceIndices(&surface_ind[0].first, &surface_ind[0].second,
                        &surface_ind[1].first, &surface_ind[1].second,
                        &surface_ind[2].first, &surface_ind[2].second);
    vector<pair<int32_t, int32_t>> nonpad_surface_ind({make_pair(-1, -1),
                                                   make_pair(-1, -1),
                                                   make_pair(-1, -1)});
    BlockNonPadSurfacesIndices(&nonpad_surface_ind[0].first,
                               &nonpad_surface_ind[0].second,
                               &nonpad_surface_ind[1].first,
                               &nonpad_surface_ind[1].second,
                               &nonpad_surface_ind[2].first,
                               &nonpad_surface_ind[2].second);
    vector<pair<int32_t, int32_t>> loop_ind({make_pair(-1, -1),
                                         make_pair(-1, -1),
                                         make_pair(-1, -1)});
    for (int i = 0; i < 3; ++i)
    {
        if (i == axis)
        {
            if (direction == 1)
                loop_ind[i].first = surface_ind[i].second;
            else
                loop_ind[i].first = surface_ind[i].first;
            loop_ind[i].second = loop_ind[i].first;
        }
        else
        {
            loop_ind[i].first = nonpad_surface_ind[i].first;
            loop_ind[i].second = nonpad_surface_ind[i].second;
        }
    }
    vector<pair<int32_t, pair<int32_t, int32_t>>> seen_labels;
    int64_t addr;
    int32_t l;
    cout << "process id " << process_id_ << " loop ind " << loop_ind[0].first << ", " << loop_ind[0].second << ", "
    << loop_ind[1].first << ", " << loop_ind[1].second << ", "
            << loop_ind[2].first << ", " << loop_ind[2].second << endl;
    if (block_ids[axis] + direction >= 0 &&
        block_ids[axis] + direction <= tissue_dims[axis] - 1 )
    {
        for (int32_t i = loop_ind[0].first; i <= loop_ind[0].second; ++i)
            for (int32_t j = loop_ind[1].first; j <= loop_ind[1].second; ++j)
                for (int32_t k = loop_ind[2].first; k <= loop_ind[2].second; ++k)
                {
                    addr = BlockVecCoorToScalar(i, j, k);
                    l = label[addr];
                    if (l > 0)
                    {
                        if (axis == 0)
                            seen_labels.push_back(make_pair(l, make_pair(j, k)));
                        else if (axis == 1)
                            seen_labels.push_back(make_pair(l, make_pair(i, k)));
                        else
                            seen_labels.push_back(make_pair(l, make_pair(i, j)));
                    }
                }
    }
    cout << "process " << process_id_ << " face label finish" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // max number of different labels in the north / south face
    cout << "block lens: " << block_lens[0] << ", " <<  block_lens[1] << ", " <<  block_lens[2] << ", " <<block_lens[axis] << endl;
    int n_max = (int)(block_lens[0] * block_lens[1] * block_lens[2] / block_lens[axis]);
    cout << "vector size " << seen_labels.size() << ", nmax = " << n_max <<endl;
    int32_t* send_buf= new int32_t[3 * n_max + 1];
    int32_t* rec_buf = new int32_t[3 * n_max + 1];
    if (!seen_labels.empty())
    {
        send_buf[0] = (int32_t)seen_labels.size();
        int i = 1;
        for (auto it = seen_labels.cbegin(); it != seen_labels.cend(); ++it)
        {
            send_buf[i] = it->first;
            send_buf[i + 1] = it->second.first;
            send_buf[i + 2] = it->second.second;
            i += 3;
        }
    }
    else
        send_buf[0] = CC_EMPTY;
    MPI_Request request;
    MPI_Status status;
    int tag = 10;
    int target_pid, src_pid;
    int n_worker = min(n_proc_, required_worker_num_);
    bool is_receiver = false;

    if (block_ids[axis] - direction >= 0 &&
        block_ids[axis] - direction <= tissue_dims[axis] - 1)
    {
        src_pid = ProcessIDFromXYZID(block_xid_ + (axis == 2 ? -direction : 0),
                                     block_yid_ + (axis == 1 ? -direction : 0),
                                     block_zid_ + (axis == 0 ? -direction : 0));
        if (src_pid >= 0 && src_pid < n_worker)
        {
            is_receiver = true;
            cout << "process " << process_id_ << " receive from process " << src_pid << endl;
            MPI_Irecv(rec_buf, 3 * n_max + 1, MPI_INT32_T, src_pid, tag, MPI_COMM_WORLD, &request);
        }
    }
    if (block_ids[axis] + direction >= 0 &&
        block_ids[axis] + direction <= tissue_dims[axis] - 1)
    {
        target_pid = ProcessIDFromXYZID(block_xid_ + (axis == 2 ? direction : 0),
                                        block_yid_ + (axis == 1 ? direction : 0),
                                        block_zid_ + (axis == 0 ? direction : 0));
        if (target_pid >= 0 && target_pid < n_worker)
        {
            cout << "process " << process_id_ << " send to process " << target_pid << endl;
            MPI_Send(send_buf, 3 * n_max + 1, MPI_INT32_T, target_pid, tag, MPI_COMM_WORLD);
        }
    }
    if (is_receiver)
    {
        MPI_Wait(&request, &status);
        ResolveNeighborFaceLabels(rec_buf, axis, direction);
    }

    delete[] send_buf;
    delete[] rec_buf;
    MPI_Barrier(MPI_COMM_WORLD);
}

void mcp3d::ParCCLabel3D::CommunicateNeighborLabels()
{
    if (!is_worker_)
        return;
    CommunicateNeighborFaceLabels(0, -1);
    CommunicateNeighborFaceLabels(0, 1);
    CommunicateNeighborFaceLabels(1, -1);
    CommunicateNeighborFaceLabels(1, 1);
    CommunicateNeighborFaceLabels(2, -1);
    CommunicateNeighborFaceLabels(2, 1);
}

void mcp3d::ParCCLabel3D::MergeGlobalLabels(int n_label_pairs,
                                            int32_t *rec_buf)
{
    if (!is_worker_)
        return;
    int other_key, other_val;
    cout << "process " << process_id_ << " in rec buf: ";
    for (int i = 0; i < 2 * n_label_pairs; i += 2)
        cout << "(, " << rec_buf[i] << ", " <<  rec_buf[i + 1] << ")    ";
    for (int i = 0; i < 2 * n_label_pairs; i += 2)
    {
        other_key = rec_buf[i];
        other_val = rec_buf[i + 1];
        if (other_val == CC_EMPTY || other_key == CC_EMPTY || other_key == other_val)
            continue;
        if (border_labels_.find(other_key) != border_labels_.end() &&
                border_labels_.find(other_val) == border_labels_.end())
        {
            border_labels_.insert(other_val);
            // if other key is leaf
            if (global_eq_labels_[other_key] != CC_EMPTY)
                global_eq_labels_[other_val] = global_eq_labels_[other_key];
            else
                global_eq_labels_[other_val] = other_key;
        }
        // both other key and other value already in border label set
        else if (border_labels_.find(other_key) != border_labels_.end() &&
                border_labels_.find(other_val) != border_labels_.end())
        {
            if (global_eq_labels_.at(other_key) == other_val || global_eq_labels_.at(other_val) == other_key)
                continue;
            int32_t par_other_key = global_eq_labels_.at(other_key);
            int32_t par_other_val = global_eq_labels_.at(other_val);
            // both head. move other key and its children to become children of other val
            if (par_other_key == par_other_val && par_other_key == CC_EMPTY)
            {
                for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
                {
                    if (it->second == other_key)
                        it->second = other_val;
                }
                global_eq_labels_[other_key] = other_val;
            }
            if (par_other_key != par_other_val)
            {
                // other key is head, other value is not. move other value and its parent and sibling to become children of other key
                if (par_other_key == CC_EMPTY)
                {
                    for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
                    {
                        if (it->second == par_other_val)
                            it->second = other_key;
                    }
                    global_eq_labels_[other_val] = other_key;
                    global_eq_labels_[par_other_val] = other_key;
                }
                else if (par_other_val == CC_EMPTY)
                {
                    for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
                    {
                        if (it->second == par_other_key)
                            it->second = other_val;
                    }
                    global_eq_labels_[other_key] = other_val;
                    global_eq_labels_[par_other_key] = other_val;
                }
                // neither is head. move other key and its parent and sibling to become children of other val's parent
                else
                {
                    for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
                    {
                        if (it->second == par_other_key)
                            it->second = par_other_val;
                    }
                    global_eq_labels_[par_other_key] = par_other_val;
                    global_eq_labels_[other_key] = par_other_val;
                }
            }
        }
        else if (border_labels_.find(other_key) == border_labels_.end() &&
                border_labels_.find(other_val) == border_labels_.end())
        {
            border_labels_.insert(other_key);
            border_labels_.insert(other_val);
            global_eq_labels_[other_key] = other_val;
            global_eq_labels_[other_val] = CC_EMPTY;
        }
        else
        {
            border_labels_.insert(other_key);
            // other val is not head, move other key to be sibling of other val
            if (global_eq_labels_[other_val] != CC_EMPTY)
                global_eq_labels_[other_key] = global_eq_labels_[other_val];
            else
                global_eq_labels_[other_key] = other_val;
        }
    }
    cout << "merged global eq labels: ";
    for (auto it = global_eq_labels_.begin(); it != global_eq_labels_.end(); ++it)
        cout << "(" << it->first << ", " << it->second << ")  ";
    cout << endl;
}

void mcp3d::ParCCLabel3D::GlobalLabelReduction()
{
    if (!is_worker_)
        return;

    int src_id, tar_id, tag = 10;
    int n_label_pairs = 0;
    int32_t* send_buf = nullptr;
    int32_t* rec_buf = nullptr;
    MPI_Status status;

    // first trim this tree to have number of leafs as power of 2
    int level_max = (int)floor(log2(min(n_proc_, required_worker_num_)));
    int delta = (int)pow(2, level_max);

    if (process_id_ >= delta)
    {
        tar_id = process_id_ - delta;
        n_label_pairs = (int)global_eq_labels_.size();
        MPI_Send(&n_label_pairs, 1, MPI_INT, tar_id, tag, MPI_COMM_WORLD);
        if (n_label_pairs > 0)
        {
            send_buf = new int32_t[2 * n_label_pairs];
            int i = 0;
            for (auto it = global_eq_labels_.cbegin();
                 it != global_eq_labels_.cend(); ++it)
            {
                send_buf[i] = it->first;
                send_buf[i + 1] = it->second;
                i += 2;
            }
            MPI_Send(send_buf, 2 * n_label_pairs, MPI_INT32_T, tar_id, tag, MPI_COMM_WORLD);
        }
    }
    else if (process_id_ + delta < min(n_proc_, required_worker_num_))
    {
        src_id = process_id_ + delta;
        MPI_Recv(&n_label_pairs, 1, MPI_INT, src_id, tag, MPI_COMM_WORLD, &status);
        if (n_label_pairs > 0)
        {
            rec_buf = new int32_t[2 * n_label_pairs];
            MPI_Recv(rec_buf, 2 * n_label_pairs, MPI_INT32_T, src_id, tag, MPI_COMM_WORLD, &status);
            MergeGlobalLabels(n_label_pairs, rec_buf);
        }
    }
    cout << "tree trimmed" << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // tree communication
    int level_base;
    for (int level = 1; level <= level_max; ++level)
    {
        delete[] rec_buf;
        rec_buf = nullptr;
        delete[] send_buf;
        send_buf = nullptr;
        level_base = (int)pow(2, level);
        if (process_id_ < delta && process_id_ % level_base == level_base / 2)
        {
            tar_id = process_id_ / level_base * level_base;
            n_label_pairs = (int)global_eq_labels_.size();
            MPI_Send(&n_label_pairs, 1, MPI_INT, tar_id, tag, MPI_COMM_WORLD);
            if (n_label_pairs > 0)
            {
                send_buf = new int32_t[2 * n_label_pairs];
                int i = 0;
                for (auto it = global_eq_labels_.cbegin();
                     it != global_eq_labels_.cend(); ++it)
                {
                    send_buf[i] = it->first;
                    send_buf[i + 1] = it->second;
                    i += 2;
                }
                MPI_Send(send_buf, 2 * n_label_pairs, MPI_INT32_T, tar_id, tag,
                         MPI_COMM_WORLD);
            }
        }
        else if (process_id_ < delta && process_id_ % level_base == 0)
        {
            MPI_Recv(&n_label_pairs, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            if (n_label_pairs > 0)
            {
                rec_buf = new int32_t[2 * n_label_pairs];
                MPI_Recv(rec_buf, 2 * n_label_pairs, MPI_INT32_T, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
                MergeGlobalLabels(n_label_pairs, rec_buf);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    delete[] rec_buf;
    rec_buf = nullptr;
    delete[] send_buf;
    send_buf = nullptr;
    // owner of full equivalence data prepare send buffer
    if (process_id_ == 0)
    {
        n_label_pairs = (int)global_eq_labels_.size();
        if (n_label_pairs > 0)
        {
            send_buf = new int32_t[2 * n_label_pairs];
            int i = 0;
            for (auto it = global_eq_labels_.cbegin(); it != global_eq_labels_.cend(); ++it)
            {
                send_buf[i] = it->first;
                send_buf[i + 1] = it->second;
                i += 2;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // broadcast number of label pairs
    MPI_Bcast(&n_label_pairs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (process_id_ != 0)
        send_buf = new int32_t[2 * n_label_pairs];
    MPI_Barrier(MPI_COMM_WORLD);

    // broadcast full equivalence data
    cout << "process id = " << process_id_ << ", n_label_pairs = " << n_label_pairs << endl;
    MPI_Bcast(send_buf, 2 * n_label_pairs, MPI_INT32_T, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (process_id_ != 0)
    {
        for (int i = 0; i < 2 * n_label_pairs; i += 2)
            global_eq_labels_[send_buf[i]] = send_buf[i + 1];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] send_buf;
    delete[] rec_buf;
}

int32_t mcp3d::ParCCLabel3D::CountGlobalCC()
{
    unordered_set<int32_t> resolved_labels (final_labels_.cbegin(), final_labels_.cend());
    int n_global_cc_max =  global_n_cc_;
    int32_t* send_buf;
    int32_t* rec_buf;
    send_buf = new int32_t[n_global_cc_max + 1];
    rec_buf = new int32_t[n_global_cc_max + 1];
    int i = 0;
    for (auto it = resolved_labels.cbegin(); it != resolved_labels.cend(); ++it)
        send_buf[i++] = *it;
    send_buf[i] = CC_EMPTY;

    int level_max = (int)floor(log2(min(n_proc_, required_worker_num_)));
    cout << "level max = " << level_max << endl;
    int delta = (int)pow(2, level_max);
    int tar_id, src_id;
    int tag = 10;
    MPI_Status status;
    // trim tree to power of two leaf
    if (process_id_ >= delta)
    {
        tar_id = process_id_ - delta;
        MPI_Send(send_buf, n_global_cc_max + 1, MPI_INT32_T, tar_id, tag, MPI_COMM_WORLD);
    }
    else if (process_id_ + delta < min(n_proc_, required_worker_num_))
    {
        src_id = process_id_ + delta;
        if (src_id < min(n_proc_, required_worker_num_))
        {
            MPI_Recv(rec_buf, n_global_cc_max + 1, MPI_INT32_T, src_id, tag, MPI_COMM_WORLD, &status);
            for (i = 0; i < n_global_cc_max + 1; ++i)
            {
                if (rec_buf[i] == CC_EMPTY)
                    break;
                resolved_labels.insert(rec_buf[i]);
            }
            i = 0;
            for (auto it = resolved_labels.cbegin(); it != resolved_labels.cend(); ++it)
                send_buf[i++] = *it;
            send_buf[i] = CC_EMPTY;
        }
    }
    cout << "process id " << process_id_ << "tree trimmed counting" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    int level_base;
    for (int level = 1; level <= level_max; ++level)
    {
        level_base = (int) pow(2, level);
        if (process_id_ < delta && process_id_ % level_base == level_base / 2)
        {
            tar_id = process_id_ / level_base * level_base;
            MPI_Send(send_buf, n_global_cc_max + 1, MPI_INT, tar_id, tag,
                     MPI_COMM_WORLD);
            cout << "process " << process_id_ << " send to " << tar_id << endl;
        }
        else if (process_id_ < delta && process_id_ % level_base == 0)
        {
            MPI_Recv(rec_buf, n_global_cc_max + 1, MPI_INT, MPI_ANY_SOURCE,
                     tag, MPI_COMM_WORLD, &status);
            for (i = 0; i < n_global_cc_max + 1; ++i)
            {
                if (rec_buf[i] == CC_EMPTY)
                    break;
                resolved_labels.insert(rec_buf[i]);
            }
            i = 0;
            for (auto it = resolved_labels.cbegin(); it != resolved_labels.cend(); ++it)
                send_buf[i++] = *it;
            send_buf[i] = CC_EMPTY;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (process_id_ == 0)
    {
        global_n_cc_ = 0;
        for (auto it = resolved_labels.cbegin(); it != resolved_labels.cend(); ++it)
            ++global_n_cc_;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&global_n_cc_, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
    return global_n_cc_;
}

void mcp3d::ParCCLabel3D::FinalizeLabel()
{
    if (!is_worker_)
        return;
    int32_t repl, l;
    for (int64_t i = 0; i < block_total_voxels_; ++i)
    {
        if (label[i] == 0)
            continue;
        l = label[i];
        if (global_eq_labels_.find(l) != global_eq_labels_.end())
        {
            repl = global_eq_labels_.at(l);
            if (repl != CC_EMPTY && global_eq_labels_[repl] != CC_EMPTY)
            {
                string msg("error in label equivalence mapping, parent l = "
                           + to_string(repl) + ", grandparent of l = "
                           + to_string(global_eq_labels_[repl]));
                MCP3D_RUNTIME_ERROR(msg)
            }
            if (repl != CC_EMPTY)
                l = repl;
            label[i] = l;
            //cout << "process " << process_id_ << ", label[i] = " << label[i] << ", global_eq_labels.at(label[i]) = " << global_eq_labels.at(label[i]) << endl;
        }
        final_labels_.insert(label[i]);
    }
    label_finalized_ = true;
}

void mcp3d::ParCCLabel3D::BuildLocalComponentRecords()
{
    if (!is_worker_)
        return;
    if (!label_finalized_) MCP3D_RUNTIME_ERROR(
            "can not build component records before global labels are finalized")
    if (final_labels_.empty())
        return;
    vector<pair<int32_t, int32_t>> nonpad_surface_ind({make_pair(-1, -1),
                                                   make_pair(-1, -1),
                                                   make_pair(-1, -1)});
    BlockNonPadSurfacesIndices(&nonpad_surface_ind[0].first,
                               &nonpad_surface_ind[0].second,
                               &nonpad_surface_ind[1].first,
                               &nonpad_surface_ind[1].second,
                               &nonpad_surface_ind[2].first,
                               &nonpad_surface_ind[2].second);
    int64_t addr;
    int32_t l;
    for (int32_t i = nonpad_surface_ind[0].first; i <= nonpad_surface_ind[0].second; ++i)
        for (int32_t j = nonpad_surface_ind[1].first; j <= nonpad_surface_ind[1].second; ++j)
            for (int32_t k = nonpad_surface_ind[2].first; k <= nonpad_surface_ind[2].second; ++k)
            {
                addr = BlockVecCoorToScalar(i, j, k);
                l = label[addr];
                if (l > 0)
                {
                    if (component_records_.find(l) == component_records_.end())
                        component_records_.emplace(make_pair(l, ComponentRecord(l)));
                    else
                        component_records_.at(l).AddPoint(
                                Position3D<int32_t>{k, j, i});
                }
            }
}

void mcp3d::ParCCLabel3D::WriteLocalComponents(string out_dir)
{
    if (!is_worker_)
        return;
    if (!label_finalized_) MCP3D_RUNTIME_ERROR(
            "can not write component records before global labels are finalized")
    BuildLocalComponentRecords();
    int32_t l;
    string cc_filepath;
    for (auto it = final_labels_.begin(); it != final_labels_.end(); ++it)
    {
        l = *it;
        // if component contained in local block, write to disk and erase record
        if (global_eq_labels_.find(l) == global_eq_labels_.end())
        {
            component_records_.at(l).Save(out_dir);
            component_records_.erase(l);
        }
    }
}

void mcp3d::ParCCLabel3D::CommunicateComponentRecords()
{
    if (!is_worker_)
        return;
    int n_pairs = min(required_worker_num_, n_proc_) / 2;
    bool even = min(required_worker_num_, n_proc_) % 2 == 0;
}

/* incomplete
void mcp3d::ParCCLabel3D::WriteLabeledImages()
{
    if (!is_worker_)
        return;

    string img_label_dir = mcp3d::JoinPath({tissue_dir_, "label"});
    if (process_id_ == 0 && !mcp3d::IsDir(img_label_dir))
        mcp3d::MakeDirectories(img_label_dir);

    vector<pair<int32_t, int32_t>> nonpad_surface_ind({make_pair(-1, -1),
                                                   make_pair(-1, -1),
                                                   make_pair(-1, -1)});
    BlockNonPadSurfacesIndices(&nonpad_surface_ind[0].first,
                               &nonpad_surface_ind[0].second,
                               &nonpad_surface_ind[1].first,
                               &nonpad_surface_ind[1].second,
                               &nonpad_surface_ind[2].first,
                               &nonpad_surface_ind[2].second);
    vector<pair<int32_t, int32_t>> surface_ind({make_pair(-1, -1),
                                            make_pair(-1, -1),
                                            make_pair(-1, -1)});
    BlockSurfaceIndices(&surface_ind[0].first,
                        &surface_ind[0].second,
                        &surface_ind[1].first,
                        &surface_ind[1].second,
                        &surface_ind[2].first,
                        &surface_ind[2].second);
    int z_delta = nonpad_surface_ind[0].first - surface_ind[0].first;
    tdata_t buf = _TIFFmalloc(sizeof(int32_t) *
                              (nonpad_surface_ind[1].second - nonpad_surface_ind[1].first + 1) *
                              (nonpad_surface_ind[2].second - nonpad_surface_ind[2].first + 1));
    int64_t addr_buf, addr_block;
    for (int32_t i = nonpad_surface_ind[0].first; i <= nonpad_surface_ind[0].second; ++i)
    {
        addr_buf = 0;
        for (int32_t j = nonpad_surface_ind[1].first; j <= nonpad_surface_ind[1].second; ++j)
            for (int32_t k = nonpad_surface_ind[2].first; k <= nonpad_surface_ind[2].second; ++k)
            {
                addr_block = BlockVecCoorToScalar(i, j, k);
                ((int32_t*)buf)[addr_buf++] = label[addr_block];
            }
    }
    _TIFFfree(buf);
} */

int main(int argc, char** argv)
{
    if (argc < 2) MCP3D_INVALID_ARGUMENT(
            "usage: par_connected_component_3d tissue_dir [block img width] [block img height]")
    string tissue_dir = argv[1];
    int img_w, img_h, foreground;
    if (argc == 2)
    {
        img_w = MPI_CC3D_W;
        img_h = MPI_CC3D_H;
        foreground = 1;
    }
    if (argc == 3)
    {
        img_w = atoi(argv[2]);
        img_h = MPI_CC3D_H;
        foreground = 1;
    }
    if (argc  == 4)
    {
        img_w = atoi(argv[2]);
        img_h = atoi(argv[3]);
        foreground = 1;
    }
    if (argc > 4)
    {
        img_w = atoi(argv[2]);
        img_h = atoi(argv[3]);
        foreground = atoi(argv[4]);
    }
    double io_time, local_label_time, global_label_resolve_time, reassign_time, count_time;
    MPI_Init(nullptr, nullptr);
    mcp3d::ParCCLabel3D par_cc_label(tissue_dir, img_w, img_h, foreground);
    auto start = chrono::high_resolution_clock::now();
    par_cc_label.ReadInput();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> elapse = end - start;
    io_time = elapse.count();
    par_cc_label.SetForeground();
    start = chrono::high_resolution_clock::now();
    par_cc_label.LabelLocal();
    end = chrono::high_resolution_clock::now();
    elapse = end - start;
    local_label_time = elapse.count();
    par_cc_label.GlobalizeLabel();
    start = chrono::high_resolution_clock::now();
    par_cc_label.CommunicateNeighborLabels();
    par_cc_label.GlobalLabelReduction();
    end = chrono::high_resolution_clock::now();
    elapse = end - start;
    global_label_resolve_time = elapse.count();
    start = chrono::high_resolution_clock::now();
    par_cc_label.FinalizeLabel();
    end = chrono::high_resolution_clock::now();
    elapse = end - start;
    reassign_time = elapse.count();
    start = chrono::high_resolution_clock::now();
    int global_cc_n = par_cc_label.CountGlobalCC();
    end = chrono::high_resolution_clock::now();
    elapse = end - start;
    count_time = elapse.count();
    string par_cc_outdir = mcp3d::JoinPath({mcp3d::parallel_module_dir(),
                                            "par_cc3d"});
    mcp3d::RemovePath(par_cc_outdir);
    mcp3d::MakeDirectories(par_cc_outdir);

    MPI_Barrier(MPI_COMM_WORLD);
    fstream of(mcp3d::JoinPath({par_cc_outdir, "par_cc3d_" + to_string(
            par_cc_label.process_id()) + ".txt"}), fstream::out);
    cout << "global nn number = " << global_cc_n << endl;
    of << "global nn number = " << global_cc_n << endl;
    cout << "number of process: " << par_cc_label.n_proc() << endl;
    of << "number of process: " << par_cc_label.n_proc() << endl;
    cout << "tissue dimension (layers, height, width) = ("
         << par_cc_label.tissue_lens()[0] << ", "
         << par_cc_label.tissue_lens()[1] << ", "
         << par_cc_label.tissue_lens()[2] << ")" << endl;
    of << "tissue dimension (layers, height, width) = ("
         << par_cc_label.tissue_lens()[0] << ", "
         << par_cc_label.tissue_lens()[1] << ", "
         << par_cc_label.tissue_lens()[2] << ")" << endl;
    cout << "block dimension (layers, height, width) = ("
         << par_cc_label.block_dims()[0] << ", "
         << par_cc_label.block_dims()[1] << ", "
         << par_cc_label.block_dims()[2] << ")" << endl;
    of << "block dimension (layers, height, width) = ("
         << par_cc_label.block_dims()[0] << ", "
         << par_cc_label.block_dims()[1] << ", "
         << par_cc_label.block_dims()[2] << ")" << endl;
    cout << "time to read input: " << io_time << " ms" << endl;
    of << "time to read input: " << io_time << " ms" << endl;
    cout << "time to label local image " << local_label_time << " ms" << endl;
    of << "time to label local image " << local_label_time << " ms" << endl;
    cout << "time to communiate and resolve global label " << global_label_resolve_time << " ms" << endl;
    of << "time to communiate and resolve global label " << global_label_resolve_time << " ms" << endl;
    cout << "time to reassign to resolved global label " << reassign_time << " ms" << endl;
    of << "time to reassign to resolved global label " << reassign_time << " ms" << endl;
    cout << "time to count number of components " << count_time << " ms" << endl;
    of << "time to count number of components " << count_time << " ms" << endl;
    of.close();
    MPI_Finalize();
}
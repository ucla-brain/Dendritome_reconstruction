//
// Created by muyezhu on 12/5/17.
//
#ifndef MCP3D_CONNECTED_COMPONENT_HPP
#define MCP3D_CONNECTED_COMPONENT_HPP

#include <iostream>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "common/mcp3d_common.hpp"
#include "image/mcp3d_image_constants.hpp"
#include "image/mcp3d_image.hpp"


const static int32_t CC_EMPTY = -1;

namespace mcp3d
{
/// image dimension is assumed to be 2d or 3d,
/// though nd support can be easily included
/// implementation based on He at al, 2017:
/// The connected-component labeling problem: A review of state-of-the-art algorithms
class EqLabelSet
{
public:
    explicit EqLabelSet(const std::vector<int64_t>& dimensions);
    int32_t n_cc() { return n_cc_; };
    int32_t max_cc_num()  { return max_cc_num_; };
    int32_t next_local_label()  { return next_local_label_; };
    int32_t rep(int32_t i)
    {
        #ifdef DEBUG
        if (i >= max_cc_num_)
            MCP3D_INVALID_ARGUMENT("cc index out of range")
        #endif
        return rep_[i];
    }
    int32_t next(int32_t i)
    {
        #ifdef DEBUG
        if (i >= max_cc_num_)
            MCP3D_INVALID_ARGUMENT("cc index out of range")
        #endif
        return next_[i];
    }
    int32_t last(int32_t i)
    {
        #ifdef DEBUG
        if (i >= max_cc_num_)
            MCP3D_INVALID_ARGUMENT("cc index out of range")
        #endif
        return last_[i];
    }
    int32_t AddNew();
    int32_t Merge(int32_t l1, int32_t l2);
    int32_t Merge(const std::vector<int32_t>& ls);
    int32_t Merge(const std::unordered_set<int32_t>& ls);
    void ConsecutiveLabel();
    void FormatPrintRep();
    void FormatPrintNext();
    void FormatPrintLast();
private:
    int64_t max_cc_num_;
    int32_t image_id_;
    int32_t next_local_label_, n_cc_;
    // rep_[i]: representative label of label i
    // next_[i]: next label value in same label set as label i
    // last_[i]: last label value in label set with representative label i
    std::vector<int32_t> rep_;
    std::vector<int32_t> next_;
    std::vector<int32_t> last_;
    std::vector<int32_t> size_;
};

void DecisionTreeCC3D(int32_t*** image, const std::vector<int64_t>& dimensions);

/// initialize zeros array
int32_t *InitLabel(const std::vector<int64_t> &dimensions);

void ExtractVectorCoordinatesFromScalar(int64_t addr_scalar,
                                        const std::vector<int64_t>& dimensions,
                                        int64_t* i, int64_t* j, int64_t* k);

void ConvertVectorCoordinatesToScalar(int64_t* addr_scalar,
                                      const std::vector<int64_t>& dimensions,
                                      const int64_t& i, const int64_t& j,
                                      const int64_t& k);

template <typename T = int32_t>
void AssignRunLabel(T* label, const std::vector<int64_t>& dimensions,
                    const int64_t& i0, const int64_t& j0, const int64_t& k0,
                    const int64_t& i1, const int64_t& j1, const int64_t& k1, int32_t l)
{
    int64_t run_start, run_end;
    mcp3d::ConvertVectorCoordinatesToScalar(&run_start, dimensions, i0, j0, k0);
    mcp3d::ConvertVectorCoordinatesToScalar(&run_end, dimensions, i1, j1, k1);
    for (int64_t i = run_start; i <= run_end; ++i)
        label[i] = l;
}

template <typename T = int32_t>
void ReassignRunLabel(T* label, const std::vector<int64_t>& dimensions,
                    const int64_t& i0, const int64_t& j0, const int64_t& k0,
                    const int64_t& i1, const int64_t& j1, const int64_t& k1, int32_t l)
{
    int64_t run_start, run_end;
    mcp3d::ConvertVectorCoordinatesToScalar(&run_start, dimensions, i0, j0, k0);
    mcp3d::ConvertVectorCoordinatesToScalar(&run_end, dimensions, i1, j1, k1);
    for (int64_t i = run_start; i <= run_end; ++i)
        if (label[i] > 0)
            label[i] = l;
}

template <typename T = int32_t>
int32_t CheckRunNeighbors(T* label, const std::vector<int64_t>& dimensions,
                        EqLabelSet* label_sets,
                       const int64_t& i0, const int64_t& j0, const int64_t& k0,
                       const int64_t& i1, const int64_t& j1, const int64_t& k1)
{
    int64_t d1 = dimensions[1],
            d2 = dimensions[2];
    std::unordered_set<int32_t> seen_label;
    int32_t l;
    int64_t i_start = std::max(i0 - 1, 0l),
            j_start = std::max(j0 - 1, 0l),
            j_end = std::min(j0 + 1, d1 - 1),
            k_start = std::max(k0 - 1, 0l),
            k_end = std::min(k1 + 1, d2 - 1);
    int64_t addr_ijk;
    if (i_start < i0)
    {
        ConvertVectorCoordinatesToScalar(&addr_ijk, dimensions, i_start, j_start, k_start);
        for (int64_t j = j_start; j <= j_end; ++j)
        {
            for (int64_t k = k_start; k <= k_end; ++k)
            {
                l = label[addr_ijk++];
                if (l > 0)
                    seen_label.insert(l);
            }
            addr_ijk += (d2 - (k_end - k_start + 1));
        }

    }
    if (j_start < j0)
    {
        ConvertVectorCoordinatesToScalar(&addr_ijk, dimensions, i0, j_start, k_start);
        for (int64_t k = k_start; k <= k_end; ++k)
        {
            l = label[addr_ijk++];
            if (l > 0)
                seen_label.insert(l);
        }
    }
    if (seen_label.empty())
        return CC_EMPTY;
    else if (seen_label.size() == 1)
        return *seen_label.begin();
    else
    {
        l = label_sets->Merge(seen_label);
        //if (i_start < i0 && j0 > j_end)
            //ReassignRunLabel(label, dimensions, i_start, j_end, k_start, i_start, j_end, k_end, l);
        return l;
    }
}

template <typename T = int32_t>
void ProcessRun(T* label, const std::vector<int64_t> &dimensions,
                EqLabelSet* label_sets,
                const int64_t& i0, const int64_t& j0, const int64_t& k0,
                const int64_t& i1, const int64_t& j1, const int64_t& k1)
{
    T neighbor_l, prov_l;
    neighbor_l = (T)(CheckRunNeighbors<T>(label, dimensions, label_sets,
                                      i0, j0, k0, i1, j1, k1));
    if (neighbor_l == CC_EMPTY)
    {
        prov_l = (T)(label_sets->AddNew());
        AssignRunLabel<T>(label, dimensions, i0, j0, k0, i1, j1, k1, prov_l);
    }
    else
    {
        prov_l = neighbor_l;
        AssignRunLabel<T>(label, dimensions, i0, j0, k0, i1, j1, k1, prov_l);
    }
}

template <typename T1 = uint8_t, typename T2 = int32_t>
void RunLabelCC3DFirstScan(T1* input, T2* label,
                           const std::vector<int64_t> &dimensions,
                           mcp3d::EqLabelSet* label_sets,
                           std::vector<int64_t>& runs)
{
    if (dimensions.size() != 3 ||
        dimensions[0] <= 0 || dimensions[1] <= 0 || dimensions[2] <= 0) MCP3D_INVALID_ARGUMENT(
            "invalid 3d image dimensions")
    bool mid_run;
    int64_t d0 = dimensions[0],
            d1 = dimensions[1],
            d2 = dimensions[2];
    int64_t n = d0 * d1 * d2, addr_ijk = 0;
    int64_t i0, j0, k0, i1, j1, k1;
    int64_t run_start, run_end;
    T1 pixel_val;
    for (int64_t addr_ij = 0; addr_ij < n - d2; addr_ij += d2)
    {
        mid_run = false;
        for (int64_t addr_k = 0; addr_k < dimensions[2]; ++addr_k)
        {
            addr_ijk = addr_ij + addr_k;
            pixel_val = input[addr_ijk];
            if (!mid_run && pixel_val > 0)
            {
                mid_run = true;
                run_start = addr_ijk;
                runs.push_back(run_start);
                ExtractVectorCoordinatesFromScalar(run_start, dimensions,
                                                   &i0, &j0, &k0);
            }
            else if (mid_run && pixel_val == 0)
            {
                mid_run = false;
                run_end = addr_ijk - 1;
                ExtractVectorCoordinatesFromScalar(run_end, dimensions,
                                                   &i1, &j1, &k1);
                runs.push_back(run_end);
                ProcessRun<T2>(label, dimensions, label_sets,
                               i0, j0, k0, i1, j1, k1);
            }
        }
        if (mid_run && pixel_val > 0)
        {
            run_end = addr_ijk;
            ExtractVectorCoordinatesFromScalar(run_end, dimensions,
                                               &i1, &j1, &k1);
            runs.push_back(run_end);
            ProcessRun<T2>(label, dimensions, label_sets,
                           i0, j0, k0, i1, j1, k1);
        }
    }
}

template <typename T = int32_t>
void RunLabelCC3DSecondScan(T * label, const std::vector<int64_t> &dimensions,
                            mcp3d::EqLabelSet* label_sets,
                            const std::vector<int64_t>& runs)
{
    if (dimensions.size() != 3 ||
        dimensions[0] <= 0 || dimensions[1] <= 0 || dimensions[2] <= 0) MCP3D_INVALID_ARGUMENT(
            "invalid 3d image dimensions")
    int64_t run_start, run_end;
    int32_t provl, l;
    for (size_t i = 0; i < runs.size(); i += 2)
    {
        run_start = runs[i];
        run_end = runs[i + 1];
        provl = label[run_start];
        l = label_sets->rep(provl);
        for (int64_t ind = run_start; ind <= run_end; ++ind)
            label[ind] = l;
    }
}

/// implementation based on He at al, 2010:
/// Two Efficient Label-Equivalence-Based Connected-Component
/// Labeling Algorithms for 3D Binary Images
/// in the paper border pixels are assumed black.
/// this assumption is NOT made here
/// the factor most relevant to runtime is the length of the z dimension
/// images should be packed in [z step level, xy plane]
/// fashion to increase performance if z step levels are limited
template <typename T1 = uint8_t, typename T2 = int32_t>
T2 RunLabelCC3D(T1* input, T2* label,
                const std::vector<int64>& dimensions, int16_t image_id = 0)
{
    if (dimensions.size() != 3 ||
       dimensions[0] <= 0 || dimensions[1] <= 0 || dimensions[2] <= 0) MCP3D_INVALID_ARGUMENT(
            "invalid 3d image dimensions")

    std::vector<int64> runs;
    mcp3d::EqLabelSet label_sets(dimensions);
    RunLabelCC3DFirstScan<T1, T2>(input, label, dimensions, &label_sets, runs);
    label_sets.ConsecutiveLabel();
    RunLabelCC3DSecondScan<T2>(label, dimensions, &label_sets, runs);
    return (T2)(label_sets.n_cc());
}

template <typename T = uint32_t>
T* InitLabel(const std::vector<int64_t> &dimensions)
{
    int64 d0 = dimensions[0],
          d1 = dimensions[1],
          d2 = dimensions[2];
    int64 n = d0 * d1 * d2;
    // allocate memory for 3D array
    T* label = new (std::nothrow) T[n];
    if (label == nullptr) MCP3D_BAD_ALLOC(
            "failed to allocate memory for label array: label_y")
    // initialize all elements to zero
    memset(label, 0, sizeof(T) * n);
    return label;
}

template<typename T1 = uint8_t, typename T2 = int32_t>
std::pair<T2, T2*> LabelCC3D(T1* input, const std::vector<int64_t>& dimensions,
                             int16 stack_id = 0, int foreground = 1,
                             const std::string& method = "run")
{
    if (!std::is_same<T1, uint8_t>::value && !std::is_same<T1, uint16_t>::value) MCP3D_DOMAIN_ERROR(
            "3d connected component labeling expect uint8 or uint16 input voxel datatype")
    if (!std::is_same<T2, int32_t>::value && !std::is_same<T2, int64_t>::value) MCP3D_DOMAIN_ERROR(
            "3d connected component labeling gives int32 or int64 output voxel datatype")
    if (dimensions.size() != 3 ||
        dimensions[0] <= 0 ||
        dimensions[1] <= 0 ||
        dimensions[2] <= 0) MCP3D_INVALID_ARGUMENT(
            "invalid 3d image dimensions")
    T2 *label = InitLabel<T2>(dimensions);
    T2 n_cc;
    if (method == "run")
        n_cc = RunLabelCC3D<T1, T2>(input, label, dimensions, stack_id);
    return std::make_pair(n_cc, label);
}

}
#endif //MCP3D_CONNECTED_COMPONENT_HPP

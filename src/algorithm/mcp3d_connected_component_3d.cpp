//
// Created by muyezhu on 12/5/17.
//
#include "algorithm/mcp3d_connected_component_3d.hpp"

using namespace std;

mcp3d::EqLabelSet::EqLabelSet(const vector<int64>& dimensions):
                              next_local_label_(1), n_cc_(0),
                              rep_({CC_EMPTY}), next_({CC_EMPTY}),
                              last_({CC_EMPTY}), size_({CC_EMPTY})
{
    if (dimensions.size() < 2 || dimensions.size() > 3) MCP3D_INVALID_ARGUMENT(
            "invalid image dimension")
    max_cc_num_ = 1;
    for (size_t i = 0; i < dimensions.size(); ++i)
        max_cc_num_ *= dimensions[i];
    if (dimensions.size() == 2)
        max_cc_num_ /= 4;
    else
        max_cc_num_ /= 27;
}

int32_t mcp3d::EqLabelSet::AddNew()
{
    // label starts at 1. new label is in its own label set
    rep_.push_back(next_local_label_);
    // next_[next_local_label_] is CC_EMPTY
    next_.push_back(CC_EMPTY);
    last_.push_back(next_local_label_);
    size_.push_back(1);
    return next_local_label_++;
}

int32_t mcp3d::EqLabelSet::Merge(int32_t l1, int32_t l2)
{
    int repl1 = rep_[l1];
    int repl2 = rep_[l2];
    if (repl1 == repl2)
        return repl1;
    int repl;
    if (size_[repl1] >= size_[repl2])
        repl = repl1;
    else
        repl = repl2;
    if (repl1 == repl)
    {
        int next = repl2;
        while (next != CC_EMPTY)
        {
            rep_[next] = repl;
            next = next_[next];
        }
        next_[last_[repl1]] = repl2;
        last_[repl1] = last_[repl2];
        size_[repl1] += size_[repl2];
        return repl1;
    }
    else
    {
        int next = repl1;
        while (next != CC_EMPTY)
        {
            rep_[next] = repl;
            next = next_[next];
        }
        next_[last_[repl2]] = repl1;
        last_[repl2] = last_[repl1];
        size_[repl2] += size_[repl1];
        return repl2;
    }
}

int32_t mcp3d::EqLabelSet::Merge(const vector<int32_t>& ls)
{
    if (ls.empty()) MCP3D_INVALID_ARGUMENT("empty equavalent label vector")
    else if (ls.size() == 1)
        return ls[0];
    else
    {
        size_t max_sizei = 0;
        int max_size = 0;
        for (size_t i = 0; i < ls.size(); ++i)
            if (size_[ls[i]] > max_size)
            {
                max_size = size_[ls[i]];
                max_sizei = i;
            }
        int32_t l = ls[max_sizei];
        for (size_t i = 0; i < ls.size(); ++i)
            Merge(l, ls[i]);
        return l;
    }
}

int32_t mcp3d::EqLabelSet::Merge(const unordered_set<int32_t>& ls)
{
    vector<int32_t> lsv = vector<int32_t>(ls.begin(), ls.end());
    return Merge(lsv);
}

void mcp3d::EqLabelSet::ConsecutiveLabel()
{
    int consec_l, next;
    for (int i = 1; i < next_local_label_; ++i)
    {
        if (rep_[i] == i)
        {
            consec_l = ++n_cc_;
            rep_[i] = consec_l;
            next = next_[i];
            while (next != CC_EMPTY)
            {
                rep_[next] = consec_l;
                next = next_[next];
            }
        }
    }
}

void mcp3d::EqLabelSet::FormatPrintRep()
{
    for (int i = 0; i < next_local_label_; ++i)
        cout << "rep[" << i << "] = " << rep_[i] << ",";
    cout << endl;
}

void mcp3d::EqLabelSet::FormatPrintNext()
{
    for (int i = 0; i < next_local_label_; ++i)
        cout << "next[" << i << "] = " << next_[i] << ",";
    cout << endl;
}

void mcp3d::EqLabelSet::FormatPrintLast()
{
    for (int i = 0; i < next_local_label_; ++i)
        cout << "last[" << i << "] = " << last_[i] << ",";
    cout << endl;
}

void mcp3d::ExtractVectorCoordinatesFromScalar(int64 addr_scalar,
                                               const vector<int64>& dimensions,
                                               int64* i, int64* j, int64* k)
{
    int64 d1 = dimensions[1], d2 = dimensions[2];
    *i = addr_scalar / (d1 * d2);
    *j = (addr_scalar / d2) % d1;
    *k = addr_scalar % d2;
}

void mcp3d::ConvertVectorCoordinatesToScalar(int64* addr_scalar,
                                             const vector<int64>& dimensions,
                                             const int64& i, const int64& j,
                                             const int64& k)
{
    int64 d1 = dimensions[1], d2= dimensions[2];
    *addr_scalar = i * d1 * d2 + j * d2 + k;
}

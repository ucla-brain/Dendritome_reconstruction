//
// Created by mzhu on 12/6/17.
//
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include "common/mcp3d_common.hpp"
#include "common/mcp3d_utility.hpp"
#include "algorithm/mcp3d_connected_component_3d.hpp"
#include "test/test_cc_label_3d_single_core.hpp"

using namespace std;

const static int TEST_LABEL_D0 = 1024;
const static int TEST_LABEL_D1 = 1024;
const static int TEST_LABEL_D2 = 200;
const static int TEST_LABEL_N = 10;

string cc_data_dir()
{
    return mcp3d::test_data_dir() + "/connected_component";
}

void mcp3d::test::NdiLabel()
{
    if (!system(nullptr)) MCP3D_RUNTIME_ERROR("no command processor available")
    string command = "python " + mcp3d::test_data_dir() + "/ndi_label.py " +
                     to_string(TEST_LABEL_N) + " " +
                     to_string(TEST_LABEL_D0) + " " +
                     to_string(TEST_LABEL_D1) + " " +
                     to_string(TEST_LABEL_D2) + " " ;
    system(command.c_str());
}

bool mcp3d::test::LabelCC3DCorrect()
{
    //mcp3d::test::NdiLabel();
    for (int i = 0; i < TEST_LABEL_N; ++i) MCP3D_ASSERT(mcp3d::IsFile(
            mcp3d::JoinPath({cc_data_dir(),
                             "random" + to_string(i) + ".npy"})))
    string npy_path, ndi_label_path;
    for (int i = 0; i < 1; ++i)
    {
       // npy_path = mcp3d::JoinPath({cc_data_dir(),
                                  //"soma" + to_string(i) + ".npy"});
        npy_path = mcp3d::JoinPath({cc_data_dir(),
                                    "random0_reshape.npy"});
        //npy_path = mcp3d::JoinPath({cc_data_dir(),
                                     //"soma0.npy"});
        ndi_label_path = mcp3d::JoinPath({cc_data_dir(),
                                          "ndi_label" + to_string(i) + ".npy"});
        //cnpy::NpyArray input = cnpy::npy_load(npy_path);
        //uint8_t* data = input.data<uint8_t>();
        //vector<int64> dimensions(input.shape.begin(), input.shape.end());
       // cout << "labeling: " << npy_path << endl;
        //auto start = chrono::high_resolution_clock::now();
        //pair<int32_t, int32_t*> label_results = mcp3d::LabelCC3D<uint8_t, int32_t>(data, dimensions, 0, 1);
        //auto end = chrono::high_resolution_clock::now();
        //chrono::duration<double, milli> elapse = end - start;
        //double duration = elapse.count();
        //cout << "finished, " << label_results.first << " components" << endl;
        //cout << "time to label: " << duration << "milliseconds" << endl;
        /*
        cnpy::NpyArray ndi_label = cnpy::npy_load(ndi_label_path);
        if (! CompareWithNDILabel<int32>(label_results.second, ndi_label.data<int32>()))
        {
            cout << "incorrect label result: " << npy_path << endl;
            return false;
        }
        else
            cout << "correct: " << npy_path << endl;*/
    }
    return true;
}

bool mcp3d::test::TestEqLabelSet()
{
    mcp3d::EqLabelSet eq(vector<int64>({30, 10, 1}));
    cout << eq.max_cc_num() << endl;
    cout << eq.n_cc() << endl;
    cout << eq.next_local_label() << endl;

    for (int i = 0; i < eq.max_cc_num()-1; ++i)
        eq.AddNew();
    eq.FormatPrintRep();
    eq.FormatPrintNext();
    eq.FormatPrintLast();
    eq.Merge(vector<int32>({1, 9}));
    cout << "eq.Merge({1, 9})" << endl;
    eq.FormatPrintRep();
    eq.FormatPrintNext();
    eq.FormatPrintLast();
    eq.Merge(vector<int32>({6,9}));
    cout << "eq.Merge({6,9})" << endl;
    eq.FormatPrintRep();
    eq.FormatPrintNext();
    eq.FormatPrintLast();
    eq.Merge(vector<int32>({2,7}));
    cout << "eq.Merge({2,7})" << endl;
    eq.FormatPrintRep();
    eq.FormatPrintNext();
    eq.FormatPrintLast();
    eq.ConsecutiveLabel();
    eq.FormatPrintRep();
    eq.FormatPrintNext();
    eq.FormatPrintLast();
    return true;
}
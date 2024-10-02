//
// Created by muyezhu on 12/10/18.
//
#include "common/mcp3d_utility.hpp"
#include "test_utility_internal.hpp"


using namespace std;

int mcp3d::test::CreateDirectoryTree(const string& starting_dir, long seed)
{
    double p_new = 0.9, p_dir = 0.1, p_hidden = 0.05;
    bool make_new, make_dir, make_hidden;
    int n_total = 1, n_max = 10000;

    default_random_engine generator(seed);
    bernoulli_distribution distribution_new(p_new);
    bernoulli_distribution distribution_dir(p_dir);
    bernoulli_distribution distribution_hidden(p_hidden);

    string current_root = starting_dir;
    string name, cmd;

    while (true)
    {
        // ensure some contents are always populated
        make_new = distribution_new(generator) || n_total == 1;
        make_dir = distribution_dir(generator);
        make_hidden = distribution_hidden(generator);
        // give an upper bound
        if (n_total == n_max)
            break;
        if (!make_new)
        {
            // termination condition
            if (current_root == starting_dir)
                break;
            current_root = mcp3d::ParentDir(current_root);
            continue;
        }
        ++n_total;
        name.clear();
        if (make_hidden)
            name.append(".");
        name.append(to_string(n_total));
        if (make_dir)
        {
            current_root.append("/");
            current_root.append(name);
            cmd = mcp3d::JoinVector<string>({"mkdir", current_root}, " ");
        }
        else
            cmd = mcp3d::JoinVector<string>({"touch", mcp3d::JoinPath({current_root, name})}, " ");
        mcp3d::SysCmdResult(cmd);
    }
    return n_total;
}

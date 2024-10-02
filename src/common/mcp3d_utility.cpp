//
// Created by mzhu on 3/13/17.
//
#include <cstdio>
#include <fstream>
#include <utility>
#include <locale>
#include <chrono>
#include <random>
#include <algorithm>
#include <iterator>
#include <regex>
#include <thread>
#include <mutex>
#include <future>
#include <boost/predef.h>
#include <boost/algorithm/string/predicate.hpp>
#include "mcp3d_utility.hpp"

using namespace std;

#if BOOST_OS_LINUX
mcp3d::DFSDirectoryWalker::DFSDirectoryWalker(const string &root_dir,
                                              bool sort_alpha):
                                                sort_alpha_(sort_alpha),
                                                tree_root_(root_dir)
{
    if (!mcp3d::IsFile(root_dir) && !mcp3d::IsDir(root_dir))
        MCP3D_RUNTIME_ERROR(root_dir + " does not exist")
    tree_stack_.emplace_back(tree_root_);
    // place root dir in fronter tree path
    tree_path_to_frontier_dir_.push_back(-1);
    // place one directory count at root
    terminal_dir_counts_.push_back(0);
    terminal_entry_counts_.push_back(1);
}

int mcp3d::DFSDirectoryWalker::TerminalDirCounts()
{
    return entries_in_terminal_dir_ == tree_stack_.back() ? terminal_dir_counts_.back() : 0;
}

int mcp3d::DFSDirectoryWalker::TerminalEntryCounts()
{
    return entries_in_terminal_dir_ == tree_stack_.back() ?
           terminal_entry_counts_.back() : CountEntries(entries_in_terminal_dir_);
}

void mcp3d::DFSDirectoryWalker::ExpandNextTerminal()
{
    // level of directories
    string frontier_dir_path;
    while (!tree_path_to_frontier_dir_.empty())
    {
        // if tree_path_to_frontier_dir_.back() == INT32_MIN,
        // indicating all files at the level has been expanded,
        // move one level up
        if (tree_path_to_frontier_dir_.back() == INT32_MIN)
        {
            tree_path_to_frontier_dir_.pop_back();
            terminal_dir_counts_.pop_back();
            terminal_entry_counts_.pop_back();
            tree_stack_.pop_back();
            entries_in_terminal_dir_.clear();
            continue;
        }

        // if tree_path_to_frontier_dir_.back() not negative (aka visited) and is not right most,
        // move one position left on same depth. encode tree path last entry as
        // -pos - 1. if tree_path_to_frontier_dir_.back() is maximum given number
        // of entries in the terminal directory, make it INT32_MIN
        // and place all contents in its parent entries_in_terminal_dir_
        // if tree_path_to_frontier_dir_.back() is negative but not
        // INT32_MIN, no expansion has occured previously, expand it
        else if (tree_path_to_frontier_dir_.back() >= 0)
        {
            if (tree_path_to_frontier_dir_.back() == terminal_entry_counts_.back() - 1)
            {
                tree_path_to_frontier_dir_.back() = INT32_MIN;
                if (tree_stack_.size() > tree_path_to_frontier_dir_.size())
                {
                    tree_stack_.pop_back();
                    terminal_dir_counts_.pop_back();
                    terminal_entry_counts_.pop_back();
                }
                entries_in_terminal_dir_ = tree_stack_.back();
                break;
            }
            else
            {
                // move right
                ++tree_path_to_frontier_dir_.back();
                // encode
                tree_path_to_frontier_dir_.back() *= -1;
                --tree_path_to_frontier_dir_.back();
            }
        }

        frontier_dir_path = ConstructTerminalDirPath();
        // expand if able
        if (mcp3d::IsDir(frontier_dir_path))
        {
            // increment level's directory count (first time seeing a directory)
            terminal_dir_counts_.back() += 1;
            // mark tree_path_to_frontier_dir_.back() to positive,
            // indicating expansion has occured
            if (tree_path_to_frontier_dir_.back() < 0)
            {
                tree_path_to_frontier_dir_.back() += 1;
                tree_path_to_frontier_dir_.back() *= -1;
            }
            // expand one depth down
            entries_in_terminal_dir_ = mcp3d::ListDirString(frontier_dir_path, sort_alpha_, true);
            if (entries_in_terminal_dir_.empty())
            // continue to next loop expansion if directory is empty
                continue;
            // stop expansion if only files are found in directory
            else if (mcp3d::IsFile(mcp3d::JoinPath({frontier_dir_path, EntryAtTreePosition((int)tree_stack_.size(), 0)})))
                break;
            // move deeper if directories found in directory. prior to expansion,
            // last entry in tree path should be -pos - 1. reconstruct later as
            // (tree_path[last entry] + 1) * -1
            else
            {
                int pos = 0;
                tree_path_to_frontier_dir_.push_back(- pos - 1);
                tree_stack_.emplace_back();
                tree_stack_.back().swap(entries_in_terminal_dir_);
                // the new level's directory will be counted when its expanded
                terminal_dir_counts_.push_back(0);
                terminal_entry_counts_.push_back(CountEntries(tree_stack_.back()));
            }
        }
        // else list same level files and leave expansion
        else
        {
            tree_path_to_frontier_dir_.back() = INT32_MIN;
            entries_in_terminal_dir_ = tree_stack_.back();
            break;
        }
    }
    // clean up exausted tree
    if (tree_path_to_frontier_dir_.size() == 1 &&
        tree_path_to_frontier_dir_.back() == INT32_MIN)
    {
        tree_stack_.clear();
        tree_path_to_frontier_dir_.clear();
        terminal_dir_counts_.clear();
        terminal_entry_counts_.clear();
        entries_in_terminal_dir_.clear();
        #ifdef VERBOSE
        MCP3D_MESSAGE("DFSDirectory walker at " + tree_root_ + " exhausted")
        #endif
    }
}

// if i = (int)tree_stack_.size() will retrieve from entries_in_terminal_dir_
string mcp3d::DFSDirectoryWalker::EntryAtTreePosition(int level, int i)
{
    if (level > (int)tree_stack_.size())
        MCP3D_RUNTIME_ERROR("tree level " + to_string(level) + " out of range")
    const string& level_files = level == -1 ? tree_stack_.back() :
                                (level == (int)tree_stack_.size() ? entries_in_terminal_dir_ : tree_stack_[level]);
    int n_newline = 0;
    size_t j = 0;
    string output;
    for( ; j < level_files.size(); ++j)
    {
        if (level_files[j] == '\n')
            ++n_newline;
        if (n_newline == i + 1)
            break;
    }
    if (j < level_files.size() - 1)
        --j;
    size_t k = j;
    while (k > 0)
        if (level_files[--k] == '\n')
            break;
    if (k > 0)
        ++k;
    return level_files.substr(k, j - k + 1);
}

string mcp3d::DFSDirectoryWalker::ConstructTerminalDirPath()
{
    string frontier_dir_path;
    for (size_t i = 0; i < tree_path_to_frontier_dir_.size(); ++i)
    {
        if (tree_path_to_frontier_dir_[i] >= 0)
        {
            if (i > 0)
                frontier_dir_path.append("/");
            frontier_dir_path.append(EntryAtTreePosition((int)i, tree_path_to_frontier_dir_[i]));
        }
        else if (tree_path_to_frontier_dir_[i] > INT32_MIN)
        {
            int real_pos = - (tree_path_to_frontier_dir_[i] + 1);
            if (i > 0)
                frontier_dir_path.append("/");
            frontier_dir_path.append(EntryAtTreePosition((int)i, real_pos));
        }
    }
    return frontier_dir_path;
}

int mcp3d::DFSDirectoryWalker::CountEntries(const string& contents)
{
    int n_entries = 0;
    for (const auto& entry_char: contents)
        if (entry_char == '\n')
            ++n_entries;
    return n_entries + 1;
}
#endif

string mcp3d::NormPathSeparatorOS(const std::string &input_path)
{
    string norm_path(input_path);
    if (BOOST_OS_WINDOWS)
        replace(norm_path.begin(), norm_path.end(), mcp3d::SeparatorChar(true), mcp3d::SeparatorChar(false));
    else
        replace(norm_path.begin(), norm_path.end(), mcp3d::SeparatorChar(false), mcp3d::SeparatorChar(true));
    return norm_path;
}

string mcp3d::NormPathDriveOS(const std::string &input_path)
{
    string norm_path(input_path);
    #if BOOST_OS_WINDOWS
        if (boost::algorithm::starts_with(norm_path, mcp3d::mcp_dir()))
           norm_path.replace(0, mcp3d::mcp_dir().size(), mcp3d::z_drive());
    #elif BOOST_OS_LINUX
        if (boost::algorithm::starts_with(norm_path, mcp3d::z_drive()))
            norm_path.replace(0, mcp3d::z_drive().size(), mcp3d::mcp_dir());
    #else
        MCP3D_RUNTIME_ERROR("only supporting Linux and Windows OS")
    #endif
    return norm_path;
}

string mcp3d::NormPathOS(const string &input_path)
{
    string norm_path(NormPathDriveOS(input_path));
    norm_path = mcp3d::NormPathSeparatorOS(input_path);
    return norm_path;
}

boost::filesystem::path mcp3d::NormBoostPathOS(const string &input_path)
{
    string norm_path(mcp3d::NormPathOS(input_path));
    return boost::filesystem::path(norm_path);
}

bool mcp3d::IsFile(const string &file_path)
{
    string norm_file_path(mcp3d::NormPathOS(file_path));
    return boost::filesystem::is_regular_file(boost::filesystem::path(norm_file_path));
}

bool mcp3d::IsDir(const string &dir_path)
{
    boost::filesystem::path p(mcp3d::NormBoostPathOS(dir_path));
    return boost::filesystem::is_directory(p);
}

bool mcp3d::EmptyDir(const std::string &dir_path)
{
    string norm_path = mcp3d::NormPathOS(dir_path);
    if (!boost::filesystem::is_directory(norm_path))
        return false;
    return boost::filesystem::is_empty(norm_path);
}

boost::filesystem::path mcp3d::ReadSymlink(const std::string &input_path)
{
    if (!boost::filesystem::is_symlink(input_path))
        return boost::filesystem::path{};
    return boost::filesystem::read_symlink(input_path);
}

boost::filesystem::path mcp3d::ResolveSymlink(const std::string &path_)
{
    string norm_path(mcp3d::NormPathOS(path_));
    vector<string> path_components = mcp3d::PathComponents(norm_path);
    boost::filesystem::path symlink_path{};
    for (const auto& path_component: path_components)
    {
        symlink_path.append(path_component);
        // follow link chain
        while (boost::filesystem::is_symlink(symlink_path))
            symlink_path = boost::filesystem::read_symlink(symlink_path);
    }
    return symlink_path;
}

pair<string, string> mcp3d::SplitBaseName(const string &file_path)
{
    return make_pair(mcp3d::ParentDir(file_path), mcp3d::Basename(file_path));
}

pair<string, string> mcp3d::SplitFirstDir(const std::string &file_path)
{
    if (file_path.empty())
        return make_pair("", "");
    string first_dir, remainder;
    first_dir = mcp3d::LeadingDir(file_path);
    remainder = file_path.substr(first_dir.size());
    remainder = mcp3d::RemoveLeadingSeparator(remainder);
    return make_pair(first_dir, remainder);
}

string mcp3d::LeadingDir(const string &file_path)
{
    if (file_path.empty())
        return "";
    string norm_file_path(mcp3d::NormPathOS(file_path));
    #if BOOST_OS_LINUX
        if (norm_file_path[0] == mcp3d::SeparatorChar())
            return "/";
        return norm_file_path.substr(0, norm_file_path.find(mcp3d::SeparatorChar()));
    #elif BOOST_OS_WINDOWS
        norm_file_path = mcp3d::RemoveLeadingSeparator(norm_file_path);
        return norm_file_path.substr(0, norm_file_path.find(mcp3d::SeparatorChar()));
    #else
        MCP3D_RUNTIME_ERROR("only supporting Linux and Windows OS")
    #endif
}

string mcp3d::ParentDir(const string &file_path)
{
    string norm_file_path(mcp3d::NormPathOS(file_path));
    norm_file_path = mcp3d::RemoveTrailingSeparator(norm_file_path);
    boost::filesystem::path d(norm_file_path);
    return d.parent_path().string();
}

string mcp3d::Basename(const string &path_)
{
    string norm_file_path(mcp3d::NormPathOS(path_));
    norm_file_path = mcp3d::RemoveTrailingSeparator(norm_file_path);
    norm_file_path = mcp3d::RemoveTrailingSeparator(norm_file_path);
    boost::filesystem::path b(norm_file_path);
    return b.filename().string();
}

vector<string> mcp3d::Basenames(const vector<string> &paths)
{
    vector<string> basenames;
    for (const auto& path: paths)
        basenames.push_back(mcp3d::Basename(path));
    return basenames;
}

string mcp3d::RemoveFileExt(const string &file_name)
{
    string s;
    size_t i = file_name.find_last_of('.');
    if (i < file_name.size())
        s = file_name.substr(0, i);
    else
        s = file_name;
    return s;
}

string mcp3d::FileExt(const string &file_name)
{
    string s;
    size_t i = file_name.find_last_of('.');
    if (i < file_name.size())
        s = file_name.substr(i);
    else
        s = "";
    return s;
}

bool mcp3d::MakeDirectories(const std::string &dir_path)
{
    if (mcp3d::IsDir(dir_path))
    {
        #ifdef VERBOSE
            MCP3D_MESSAGE(dir_path + " already exists")
        #endif
        return true;
    }
    boost::filesystem::path directory(dir_path);
    try
    {
        MCP3D_ASSERT(boost::filesystem::create_directories(directory));
        return true;
    }
    // guarding against multi-threading execution where another thread starts
    // creation after current thread existence checking but finishes earlier
    catch (...)
    {
        return mcp3d::IsDir(dir_path);
    }
}

bool mcp3d::RemovePath(const std::string& path_, bool remove_link_target)
{
    string norm_path = mcp3d::NormPathOS(path_);
    if (!boost::filesystem::exists(norm_path))
    {
        #ifdef VERBOSE
            MCP3D_MESSAGE(path_ + " does not exist. do nothing.")
        #endif
        return true;
    }
    boost::filesystem::path p(norm_path);
    bool is_simlink = boost::filesystem::is_symlink(p);
    // guarding multi-thread execution, when thread A deleted path_ prior to
    // thread B
    try
    {
        // if remove_link_target, recursively follow link chain
        if (is_simlink && remove_link_target)
        {
            boost::filesystem::path p_target(boost::filesystem::read_symlink(p));
            MCP3D_ASSERT(RemovePath(p_target.string(), true))
        }
        if (boost::filesystem::is_directory(p))
        {
            // p is symlink to directory
            if (is_simlink)
                MCP3D_ASSERT(boost::filesystem::remove(p))
            // p is directory
            else
                MCP3D_ASSERT(boost::filesystem::remove_all(p))
        }
        else
            MCP3D_ASSERT(boost::filesystem::remove(p))
        return true;
    }
    catch (...)
    {
        // check if p itself exists
        if (boost::filesystem::exists(p))
            return false;
        if (remove_link_target)
        {
            boost::filesystem::path target(mcp3d::ReadSymlink(p.string()));
            while (!target.empty())
            {
                if (boost::filesystem::exists(target))
                    return false;
                target = mcp3d::ReadSymlink(target.string());
            }
        }
        return true;
    }
}

void mcp3d::AsyncRemovePath(const std::string &path_, const string &mode)
{
    #if BOOST_OS_WINDOWS
        mcp3d::RemovePath(path_);
    #endif
    #if BOOST_OS_LINUX
    // if path_ is a file, delete it and return
    if (mcp3d::IsFile(path_))
    {
        mcp3d::RemovePath(path_);
        return;
    }
    else if (!mcp3d::IsDir(path_))
    {
        MCP3D_MESSAGE(path_ + " does not exist, do nothing")
        return;
    }

    mcp3d::DFSDirectoryWalker walker(path_, false);

    string cmd_files, cmd_dirs;
    bool exausted = false;
    int n_dirs, n_dir_chars, n_files;
    future<bool> file_complete, dir_complete;

    while (!exausted)
    {
        // non blocking io requests
        if (!cmd_files.empty())
        {
            file_complete = async(launch::async, [&](){
                SysCmdResult(cmd_files);
                return true;
            });
        }

        if (!cmd_dirs.empty())
        {
            dir_complete = async(launch::async, [&](){
                SysCmdResult(cmd_dirs);
                return true;
            });
        }


        string file_names, dir_names;
        walker.ExpandNextTerminal();

        if (walker.Exausted())
        {
            exausted = true;
            if (!cmd_files.empty())
                file_complete.get();
            if (!cmd_dirs.empty())
                dir_complete.get();
            break;
        }

        if (walker.TerminalFiles().empty())
        {
            if (!cmd_files.empty())
                file_complete.get();
            if (!cmd_dirs.empty())
                dir_complete.get();
            cmd_files.clear();
            cmd_dirs.clear();
            continue;
        }

        n_dirs = walker.TerminalDirCounts();
        n_files = walker.TerminalEntryCounts() - n_dirs;
        n_dir_chars = 0;

        if (n_dirs > 0)
        {
            int n_encountered_dirs = 0;
            for (const auto& terminal_files_char: walker.TerminalFiles())
            {
                ++n_dir_chars;
                if (terminal_files_char == '\n')
                    ++n_encountered_dirs;
                if (n_encountered_dirs == n_dirs)
                    break;
            }
        }

        if (n_files > 0)
        {
            file_names = walker.TerminalFiles().substr((size_t)n_dir_chars);
            for (auto& file_char:file_names)
                if (file_char == '\n')
                    file_char = ' ';
        }
        if (n_dirs > 0)
        {
            dir_names = walker.TerminalFiles().substr(0, (size_t)n_dir_chars);
            for (auto& dir_char:dir_names)
                if (dir_char == '\n')
                    dir_char = ' ';
        }

        if (mode == "verbose")
            cout << "removing " << file_names << endl;

        if (!cmd_files.empty())
            file_complete.get();
        if (!cmd_dirs.empty())
            dir_complete.get();

        cmd_files.clear();
        cmd_dirs.clear();
        if (n_files > 0)
        {
            cmd_files.append(mcp3d::JoinVector<string>({"cd", walker.TerminalDirPath(), "&&"}, " "));
            cmd_files.append(mcp3d::JoinVector<string>({"", "rm", file_names}, " "));
        }
        if (n_dirs > 0)
        {
            cmd_dirs.append(mcp3d::JoinVector<string>({"cd", walker.TerminalDirPath(), "&&"}, " "));
            cmd_dirs.append(mcp3d::JoinVector<string>({"", "rmdir", dir_names}, " "));
        }
    }
    boost::filesystem::remove(path_);
    #endif
}

string mcp3d::RemoveTrailingSeparator(const string &path_)
{
    string path(path_);
    while (!path.empty() && path != string(1, mcp3d::SeparatorChar()) &&
           path.back() == mcp3d::SeparatorChar() )
        path.pop_back();
    return path;
}

string mcp3d::RemoveLeadingSeparator(const string &path_)
{
    string path(path_);
    while (!path.empty() && path.front() == mcp3d::SeparatorChar())
        path = path.substr(1);
    return path;
}

string mcp3d::CollapseLeadingSeparator(const std::string &path_)
{
    if (path_[0] != mcp3d::SeparatorChar())
        return path_;
    size_t pos = 1;
    while (path_[pos] == mcp3d::SeparatorChar())
        ++pos;
    return path_.substr(pos - 1);
}

string mcp3d::DiscardLeadingComponent(const string &path_, const string &component_)
{
    if (component_.empty())
        return path_;
    string path = mcp3d::CollapseLeadingSeparator(path_);
    string component = mcp3d::CollapseLeadingSeparator(component_);
    if (path.find(component) != 0)
        return path_;
    string path_cut(path.substr(component_.size()));
    path_cut = mcp3d::RemoveLeadingSeparator(path_cut);
    path_cut = mcp3d::RemoveTrailingSeparator(path_cut);
    return path_cut;
}

vector<string> mcp3d::PathComponents(const string &path)
{
    vector<string> components;
    string input_path(path);
    while (true)
    {
        pair<string, string> splits = mcp3d::SplitFirstDir(input_path);
        string first_dir(splits.first);
        string remainder(splits.second);
        components.push_back(move(first_dir));
        if (remainder.empty())
            break;
        input_path = remainder;
    }
    return components;
}

template <>
std::string mcp3d::JoinVector<std::string>(const std::vector<std::string>& v,
                                           const std::string& d, bool brackets)
{
    std::string delim(d), output = brackets ? "[" : "";
    if (delim.empty())
        delim = " ";
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (i == 0)
            output.append(v[i]);
        else
        {
            output.append(delim);
            output.append(v[i]);
        }
    }
    if (brackets)
        output.append("]");
    return output;
}

pair<bool, string> mcp3d::NearestCommonDir(const vector<string> &paths)
{
    vector<string> common_dirs, file_paths(paths);
    string common_component;
    bool has_common_dir = true;
    while (has_common_dir)
    {
        for (size_t i = 0; i < file_paths.size(); ++i)
        {
            pair<string, string> splits = mcp3d::SplitFirstDir(file_paths[i]);
            if (splits.first.empty())
            {
                has_common_dir = false;
                break;
            }
            if (i == 0)
                common_component = splits.first;
            if (splits.first != common_component)
            {
                has_common_dir = false;
                break;
            }
            file_paths[i] = splits.second;
        }
        if (has_common_dir)
            common_dirs.push_back(common_component);
    }
    if (common_dirs.empty())
        return make_pair(false, string());
    else
        return make_pair(true, mcp3d::JoinPath(common_dirs));
}

string mcp3d::JoinPath(const vector<string> &components_)
{
    vector<string> components(components_);
    string out_path;
    for (int i = 0; i < (int)components.size(); ++i)
    {
        if (components[i].empty())
            continue;
        if (out_path.empty())
            out_path.append(mcp3d::CollapseLeadingSeparator(components[i]));
        else
        {
            out_path = mcp3d::RemoveTrailingSeparator(out_path);
            components[i] = mcp3d::RemoveLeadingSeparator(components[i]);
            if (out_path != string(1, mcp3d::SeparatorChar()))
                out_path.append(string(1, mcp3d::SeparatorChar()));
            out_path.append(components[i]);
        }
    }
    return out_path;
}

vector<string> mcp3d::FilesInDir(const string &directory, bool full_path, bool sort_files,
                                 const vector<string> &incl, int n, const vector<string> &excl)
{
    vector<string> files;
    if (!mcp3d::IsDir(directory))
        MCP3D_OS_ERROR(directory + " is not a directory");

    if (n == 0)
        return vector<string> {};

    boost::filesystem::path p(directory);
    boost::filesystem::directory_iterator it_end;
    vector<string> results;

    for(auto it = boost::filesystem::directory_iterator(p); it != it_end; ++it)
    {
        if (boost::filesystem::is_regular_file(it->status()))
        {
            string file_name = it->path().filename().string();
            bool has_incl = incl.empty(), has_excl = false;
            for (const auto& incl_str: incl)
            {
                if (file_name.find(incl_str) != string::npos)
                    has_incl = true;
            }
            for (const auto& excl_str: excl)
            {
                if (file_name.find(excl_str) != string::npos)
                    has_excl = true;
            }
            if (has_incl && !has_excl)
                results.push_back(file_name);
        }
    }
    if (sort_files)
        sort(results.begin(), results.end());
    if (full_path)
        for(auto& result: results)
            result = mcp3d::JoinPath(directory, result);
    if (n > 0 && (size_t)n <= results.size())
        return vector<string>(results.begin(), results.begin() + n);
    return results;
}

vector<string> mcp3d::FilesEndsWithInDir(const std::string &directory, const std::string &end_str, bool full_path, bool sort_files, int n)
{
    if (!mcp3d::IsDir(directory))
        MCP3D_RUNTIME_ERROR(directory + " is not a directory")
    if (n == 0)
        return vector<string> {};
    boost::filesystem::path p(directory);
    boost::filesystem::directory_iterator it_end;
    vector<string> results;
    for(auto it = boost::filesystem::directory_iterator(p); it != it_end; ++it)
    {
        if (boost::filesystem::is_regular_file(it->status()))
        {
            string file_name = it->path().filename().string();
            if (boost::algorithm::ends_with(file_name, end_str))
                results.push_back(file_name);
        }
    }
    if (sort_files)
        sort(results.begin(), results.end());
    if (full_path)
        for(auto& result: results)
            result = mcp3d::JoinPath(directory, result);
    if (n > 0 && (size_t)n <= results.size())
        return vector<string>(results.begin(), results.begin() + n);
    return results;
}

vector<string> mcp3d::DirsInDir(const string &directory_, bool full_path, bool sort_files,
                                const vector<string> &incl, int n, const vector<string> &excl)
{
    string directory(directory_);
    boost::filesystem::path dir_path(directory);
    if (!boost::filesystem::is_directory(dir_path))
        MCP3D_OS_ERROR(directory + " is not a valid directory")
    vector<string> dir_names;
    boost::filesystem::directory_iterator dir_iter(dir_path);
    string current_dir_name;
    bool fit_criteria;
    for (auto it = boost::filesystem::begin(dir_iter); it != boost::filesystem::end(dir_iter); ++it)
    {
        if (boost::filesystem::is_directory(it->path()))
        {
            fit_criteria = true;
            current_dir_name = it->path().filename().string();
            for (const string& in: incl)
                if (current_dir_name.find(in) >= current_dir_name.size())
                    fit_criteria = false;
            for (const string& ex: excl)
                if (current_dir_name.find(ex) < current_dir_name.size())
                    fit_criteria = false;
            if (fit_criteria)
                dir_names.push_back(current_dir_name);
        }
    }
    if (sort_files)
        sort(dir_names.begin(), dir_names.end());
    while (n > 0 && n < (int)dir_names.size())
        dir_names.pop_back();
    if (full_path)
        for (string& dir_name: dir_names)
            dir_name = mcp3d::JoinPath(directory, dir_name);
    return dir_names;
}

string mcp3d::PadNumStr(int64_t n, int width)
{
    MCP3D_ASSERT(n >= 0 && width > 0)
    MCP3D_ASSERT(to_string(n).size() <= (size_t)width)
    string n_string = to_string(n), zeros;
    for (size_t i = to_string(n).size(); i < (size_t)width; ++i)
        zeros.append("0");
    n_string= zeros.append(n_string);
    return n_string;
}

void mcp3d::RenameFilesToSortable(vector<string> &file_list, const string &prenumber, const string &postnumber, bool sort_files)
{
    int max = 0;
    unsigned long num_len = 0;
    string directory = mcp3d::ParentDir(file_list.at(0));
    bool prepend_dir = (prenumber.find(directory) >= prenumber.size());
    string prenumber_copy;
    if (prepend_dir)
        prenumber_copy = directory + prenumber;
    // find max number in the numbering
    for (const string& s: file_list)
    {
        string sc(s);
        size_t p1 = s.find_first_of(prenumber_copy);
        if (p1 >= sc.size())
            continue;
        sc.replace(p1, prenumber_copy.size(), "");
        // search for "" returns npos
        size_t p2 = sc.find_first_of(postnumber);
        if (p2 >= s.size())
            continue;
        sc.replace(p2, postnumber.size(), "");
        int num = stoi(sc);
        if (num > max)
            max = num;
    }
    num_len = to_string(max).size();
    // pad leading 0s to numbers if necessary
    for (uint32_t i = 0; i < file_list.size(); ++i)
    {
        string sc(file_list[i]);
        size_t p1 = file_list[i].find_first_of(prenumber_copy);
        if (p1 >= sc.size())
            continue;
        sc.replace(p1, prenumber_copy.size(), "");
        size_t p2 = sc.find_first_of(postnumber);
        if (p2 >= file_list[i].size())
            continue;
        sc.replace(p2, postnumber.size(), "");
        if (sc.size() < num_len)
        {
            string pad(num_len - sc.size(), '0');
            sc = pad.append(sc);
            string newname = prenumber_copy.append(sc).append(postnumber);
            int success = rename(file_list[i].c_str(), newname.c_str());
            if (success == 0)
                file_list[i] = newname;
            else
                cout << "unable to rename " + file_list[i] + " to " + newname + ", error: " << errno << endl;
        }
    }
    // sort of sort_files
    if (sort_files)
        sort(file_list.begin(), file_list.end());
}

vector<string> SplitStringWhitespaces(const string &in_string)
{
    vector<string> splitted;
    string holder;
    for (char c: in_string)
    {
        if (isspace(c))
        {
            if (!holder.empty())
            {
                splitted.push_back(holder);
                holder = string();
            }
        }
        else
            holder += c;
    }
    if (!holder.empty())
        splitted.push_back(holder);
    return splitted;
}

vector<string> SplitStringChar(const string &in_string, char delim)
{
    vector<string> splitted;
    string holder;
    for (char c: in_string)
    {
        if (c == delim)
        {
            splitted.push_back(holder);
            holder = string();
        }
        else
            holder += c;
    }
    if (!holder.empty())
        splitted.push_back(holder);
    return splitted;
}

vector<string> mcp3d::SplitString(const string &in_string, const string &delim)
{
    if (delim.empty())
        return SplitStringWhitespaces(in_string);
    if (delim.size() == 1)
    {
        if (isspace(delim[0]))
            return SplitStringWhitespaces(in_string);
        return SplitStringChar(in_string, delim[0]);
    }
    vector<string> splitted;
    string holder;
    size_t pos = 0, found;
    while (pos < in_string.size())
    {
        found = in_string.find(delim, pos);
        holder = found == string::npos ?
                 in_string.substr(pos) : in_string.substr(pos, found - pos);
        if (!holder.empty())
            splitted.push_back(holder);
        if (found == string::npos)
            break;
        pos = found + delim.size();
    }
    return splitted;
}

string mcp3d::StringLower(const string &input)
{
    locale loc;
    string lower_input;
    for(const char &c: input)
        lower_input += tolower(c, loc);
    return lower_input;
}

string& mcp3d::StringLower(string &input)
{
    locale loc;
    for(char &c: input)
        c = tolower(c, loc);
    return input;
}

string& mcp3d::StringLower(string &&input)
{
    locale loc;
    for(char &c: input)
        c = tolower(c, loc);
    return input;
}

string mcp3d::StripLeading(const std::string &input)
{
    string trim_leading;
    for (size_t i = 0; i < input.size(); ++i)
    {
        if (input[i] == ' ' || input[i] == '\n')
            continue;
        else
        {
            trim_leading = input.substr(i);
            break;
        }
    }
    return trim_leading;
}

string mcp3d::StripTrailing(const std::string &input)
{
    string trim_trailing;
    for (size_t i = input.size() - 1; i >= 0; --i)
    {
        if (input[i] == ' ' || input[i] == '\n')
            continue;
        else
        {
            trim_trailing = input.substr(0, i + 1);
            break;
        }
    }
    return trim_trailing;
}

string mcp3d::Strip(const std::string &input)
{
    string trim_leading = mcp3d::StripLeading(input);
    return mcp3d::StripTrailing(trim_leading);
}

#if BOOST_OS_LINUX

double ProcRAMMB(const string &mode, const string& u)
{
    string unit = mcp3d::StringLower(u);
    if (unit != "mb" && unit != "gb")
        MCP3D_INVALID_ARGUMENT("memory unit not understood. use either MB or GB")
    string search = mode == "current" ? "VmRSS" : "VmHWM";
    int process_id = ::getpid();
    cout << process_id << endl;
    string process_meminfo_path = mcp3d::JoinPath({"/proc",
                                                   to_string(process_id),
                                                   "status"});
    if (!mcp3d::IsFile(process_meminfo_path))
    {
        printf("can not find /proc/%d/status, return -1\n", process_id);
        return -1;
    }
    ifstream f;
    f.open(process_meminfo_path);
    if (!f.good())
    {
        printf("can not find /proc/%d/status, return -1\n", process_id);
        return -1;
    }
    string line, num_string;
    while (getline(f, line))
    {
        if (line.substr(0, search.size()) == search)
        {
            for (char c: line)
                if (isdigit(c))
                    num_string += c;
        }
    }
    double nkb = strtol(num_string.c_str(), 0, 10);
    f.close();
    double usage = nkb / 1024;
    if (unit == "mb")
        return usage;
    else
        return usage / 1024;
}

double mcp3d::ProcCurrRAM(const string& unit)
{
    string mode("current");
    return ProcRAMMB(mode, unit);
}

double mcp3d::ProcPeakRAM(const string& unit)
{
    string mode("peak");
    return ProcRAMMB(mode, unit);
}

vector<string> mcp3d::ListDir(const string &dir_path, bool sort_alpha,
                              bool list_directories_first, bool revert)
{
    if (list_directories_first)
        sort_alpha = true;
    string sort_flag = sort_alpha ? "-A1" : "-AU1";
    string group_first_flag = "--group-directories-first";
    string ls_dir_cmd = list_directories_first?
                        mcp3d::JoinVector<string>({"ls", sort_flag, group_first_flag, dir_path}, " ", false) :
                        mcp3d::JoinVector<string>({"ls", sort_flag, dir_path}, " ", false);

    // . and .. is always the first two entries
    string dir_contents = mcp3d::SysCmdResult(ls_dir_cmd.c_str());
    vector<string> results = mcp3d::SplitString(dir_contents, "\n");
    if (revert)
        reverse(results.begin(), results.end());
    return results;
}

string mcp3d::ListDirString(const string &dir_path, bool sort_alpha,
                            bool list_directories_first)
{
    if (list_directories_first)
        sort_alpha = true;
    string sort_flag = sort_alpha ? "-A1" : "-AU1";
    string group_first_flag = "--group-directories-first";
    string ls_dir_cmd = list_directories_first?
                        mcp3d::JoinVector<string>({"ls", sort_flag, group_first_flag, dir_path}, " ", false) :
                        mcp3d::JoinVector<string>({"ls", sort_flag, dir_path}, " ", false);
    return mcp3d::SysCmdResult(ls_dir_cmd.c_str());
}

string mcp3d::SysCmdResult(const char * const cmd, const string& mode)
{
    string result;
    std::array<char, 128> buffer;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), &pclose);
    if (!pipe)
    MCP3D_RUNTIME_ERROR("popen() for " + string(cmd) + " failed");
    if (mode == "verbose")
        MCP3D_MESSAGE("run command: " + string(cmd))
    while (fgets(buffer.data(), sizeof buffer, pipe.get()) != nullptr)
        result += buffer.data();
    if (result.back() == '\n')
        result.pop_back();
    return result;
}

#endif

double mcp3d::MemorySize(double nbytes, const std::string &unit)
{
    string u = mcp3d::StringLower(unit);
    if (u == "bytes")
        return nbytes;
    nbytes /= 1024;
    if (u == "kb")
        return nbytes;
    nbytes /= 1024;
    if (u == "mb")
        return nbytes;
    nbytes /= 1024;
    if (u == "gb")
        return nbytes;
    cout << __FILE__ << ", " << __FUNCTION__ << ": "
         << "memory unit " << unit << " not understood. returning memory size in GB" << endl;
    return nbytes;
}

int mcp3d::DefaultNumThreads()
{
    if (mcp3d::HostOnCluster())
        return 1;
    return thread::hardware_concurrency();
}

int mcp3d::ParseInt(const string& int_str)
{
    try
    {
        int value = stoi(int_str);
        return value;
    }
    catch(...)
    {
        MCP3D_INVALID_ARGUMENT("can not convert " + int_str + " to integer")
    }
}

double mcp3d::ParseDouble(const string& double_str)
{
    try
    {
        double value = stod(double_str);
        return value;
    }
    catch(...)
    {
        MCP3D_INVALID_ARGUMENT("can not convert " + double_str + " to double")
    }
}

string mcp3d::RandChars(int l)
{
    char random_chars[l];
    uniform_int_distribution<int> distribution(65, 90);
    int64_t seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    for (int i = 0; i < l; ++i)
        random_chars[i] = (char)distribution(generator);
    random_chars[l] = '\0';
    return string(random_chars);
}

string mcp3d::HostName()
{
    char host[255];
    gethostname(host, 255);
    return string(host);
}

bool mcp3d::HostOnCluster()
{
    string host_name(mcp3d::HostName());
    regex cluster_host_pattern("c2[0-9]{3,3}");
    smatch m;
    regex_match(host_name, m, cluster_host_pattern);
    return !m.empty();
}

#if MCP3D_MPI_BUILD

void mcp3d::AsyncRemovePathMPI(const string &path_, MPI_Comm comm,
                               const string &mode)
{
    // if path_ is a file, delete it and return
    if (mcp3d::IsFile(path_))
    {
        mcp3d::RemovePath(path_);
        return;
    }
    else if (!mcp3d::IsDir(path_))
    {
        MCP3D_MESSAGE(path_ + " does not exist, do nothing")
        return;
    }

    mcp3d::DFSDirectoryWalker walker(path_, false);
    MPI_Barrier(comm);

    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (size == 1)
        return mcp3d::AsyncRemovePath(path_);

    string cmd_files, cmd_dirs;
    bool exausted = false;
    int n_dirs, n_files;
    future<bool> file_complete, dir_complete;

    int n_entries_per_process = 0;
    string file_path;
    while (!exausted)
    {
        // non blocking io requests
        if (rank > 0 && !cmd_files.empty())
        {
            file_complete = async(launch::async, [&](){
                SysCmdResult(cmd_files);
                return true;
            });
        }

        if (rank == 0 && !cmd_dirs.empty())
        {
            dir_complete = async(launch::async, [&](){
                SysCmdResult(cmd_dirs);
                return true;
            });
        }

        string file_names, dir_names;

        walker.ExpandNextTerminal();
        MPI_Barrier(comm);

        if (walker.Exausted())
        {
            exausted = true;
            if (rank > 0 && !cmd_files.empty())
                file_complete.get();
            if (rank == 0 && !cmd_dirs.empty())
                dir_complete.get();
            MPI_Barrier(comm);
            break;
        }

        if (walker.TerminalFiles().empty())
        {
            if (rank > 0 && !cmd_files.empty())
                file_complete.get();
            if (rank == 0 && !cmd_dirs.empty())
                dir_complete.get();
            MPI_Barrier(comm);
            cmd_files.clear();
            cmd_dirs.clear();
            continue;
        }

        n_files = walker.TerminalEntryCounts() - walker.TerminalDirCounts();
        n_dirs = walker.TerminalDirCounts();
        // rank 0 only do rmdir
        if (rank == 0)
        {
            if (n_dirs > 0)
            {
                const string& terminal_files = walker.TerminalFiles();
                int n_dirs_seen = 0;
                size_t  i = 0;
                for (const auto& dir_char: terminal_files)
                {
                    if (dir_char == '\n')
                        ++n_dirs_seen;
                    ++i;
                    if (n_dirs_seen == n_dirs)
                        break;
                }
                dir_names = terminal_files.substr(0, i);
                for (auto& dir_char: dir_names)
                    if (dir_char == '\n')
                        dir_char = ' ';
            }
        }
        else
        {
            n_entries_per_process = n_files / (size - 1) + (int)(n_files % (size - 1) != 0);
            const string& terminal_files = walker.TerminalFiles();
            int start = n_entries_per_process * (rank - 1) + n_dirs;
            if (n_files > 0 && start < walker.TerminalEntryCounts())
            {
                int n_entries_seen = 0;
                // first character at i, last at j
                size_t i = 0, j = 0;
                for ( ; i < terminal_files.size(); ++i)
                {
                    if (terminal_files[i] == '\n')
                        ++n_entries_seen;
                    if (n_entries_seen == start)
                        break;
                }
                // single file no directory case

                j = i + 1;
                n_entries_seen = 0;
                for ( ; j < terminal_files.size(); ++j)
                {
                    if (terminal_files[j] == '\n')
                        ++n_entries_seen;
                    if (n_entries_seen == n_entries_per_process)
                        break;
                }

                if (i > 0)
                    ++i;
                if (j < terminal_files.size())
                    --j;
                file_names = terminal_files.substr(i, j - i + 1);

                for (auto& file_char: file_names)
                    if (file_char == '\n')
                        file_char = ' ';
            }
        }

        file_path = walker.TerminalDirPath();

        if (mode == "verbose")
            cout << "removing " << file_names << endl;

        if (rank > 0 && !cmd_files.empty())
            file_complete.get();
        if (rank == 0 && !cmd_dirs.empty())
            dir_complete.get();

        MPI_Barrier(comm);

        cmd_files.clear();
        cmd_dirs.clear();
        if (rank > 0 && !file_names.empty())
        {
            cmd_files.append(mcp3d::JoinVector<string>({"cd", walker.TerminalDirPath(), "&&"}, " "));
            cmd_files.append(mcp3d::JoinVector<string>({"", "rm", file_names}, " "));
        }
        if (rank == 0 && !dir_names.empty())
        {
            cmd_dirs.append(mcp3d::JoinVector<string>({"cd", walker.TerminalDirPath(), "&&"}, " "));
            cmd_dirs.append(mcp3d::JoinVector<string>({"", "rmdir", dir_names}, " "));
        }
    }
    try
    {
        boost::filesystem::remove(path_);
    }
    catch (...)
    {
        if (mcp3d::IsDir(path_))
            MCP3D_RETHROW(current_exception())
    }
}

#endif //MCP3D_MPI_BUILD





//
// Created by muyezhu on 3/1/18.
//
#include <fstream>
#include <algorithm>
#include <chrono>
#include <random>
#include <thread>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include "common/mcp3d_paths.hpp"
#include "common/mcp3d_utility.hpp"
#include "test_utility_internal.hpp"

using namespace std;


TEST(Utility, RemoveTrailingSeparator)
{
    EXPECT_EQ("", mcp3d::RemoveTrailingSeparator(""));
    EXPECT_EQ("/", mcp3d::RemoveTrailingSeparator("///"));
    EXPECT_EQ("/", mcp3d::RemoveTrailingSeparator("/"));
    EXPECT_EQ(".", mcp3d::RemoveTrailingSeparator(".///"));
    EXPECT_EQ("foo", mcp3d::RemoveTrailingSeparator("foo///"));
    EXPECT_EQ("foo", mcp3d::RemoveTrailingSeparator("foo/"));
}

TEST(Utility, RemoveLeadingSeparator)
{
    EXPECT_EQ("", mcp3d::RemoveLeadingSeparator(""));
    EXPECT_EQ(".///", mcp3d::RemoveLeadingSeparator(".///"));
    EXPECT_EQ("", mcp3d::RemoveLeadingSeparator("////"));
    EXPECT_EQ("foo", mcp3d::RemoveLeadingSeparator("//foo"));
}

TEST(Utility, CollapseLeadingSeparator)
{
    EXPECT_EQ("", mcp3d::CollapseLeadingSeparator(""));
    EXPECT_EQ("/", mcp3d::CollapseLeadingSeparator("/"));
    EXPECT_EQ("/", mcp3d::CollapseLeadingSeparator("///"));
    EXPECT_EQ("/foo", mcp3d::CollapseLeadingSeparator("/foo"));
    EXPECT_EQ("/foo", mcp3d::CollapseLeadingSeparator("///foo"));
}

TEST(Utility, ParentDir)
{
    EXPECT_EQ("", mcp3d::ParentDir(""));
    EXPECT_EQ(".", mcp3d::ParentDir("./foo"));
    EXPECT_EQ(".", mcp3d::ParentDir("./foo/"));
    EXPECT_EQ(".", mcp3d::ParentDir(".//foo//"));
    EXPECT_EQ("/", mcp3d::ParentDir("/foo"));
    EXPECT_EQ("/", mcp3d::ParentDir("/foo/"));
    EXPECT_EQ("", mcp3d::ParentDir("foo"));
    EXPECT_EQ("", mcp3d::ParentDir("./"));
    EXPECT_EQ("", mcp3d::ParentDir("/"));
    EXPECT_EQ("", mcp3d::ParentDir("///"));
    EXPECT_EQ("/foo", mcp3d::ParentDir("/foo/bar"));
    EXPECT_EQ("/foo", mcp3d::ParentDir("/foo/bar/"));
    EXPECT_EQ("foo", mcp3d::ParentDir("foo/bar"));
    EXPECT_EQ("foo", mcp3d::ParentDir("foo/bar/"));
}

TEST(Utility, BaseName)
{
    EXPECT_EQ("", mcp3d::ParentDir(""));
    EXPECT_EQ("bar", mcp3d::Basename("/foo/bar/"));
    EXPECT_EQ("bar", mcp3d::Basename("/foo/bar"));
    EXPECT_EQ("foo", mcp3d::Basename("/foo/"));
    EXPECT_EQ("foo", mcp3d::Basename("/foo"));
    EXPECT_EQ("foo", mcp3d::Basename("foo"));
    EXPECT_EQ("/", mcp3d::Basename("/"));
    EXPECT_EQ("/", mcp3d::Basename("///"));
    EXPECT_EQ(".", mcp3d::Basename("./"));
}

TEST(Utility, LeadingDir)
{
    EXPECT_EQ("", mcp3d::LeadingDir(""));
    EXPECT_EQ("/", mcp3d::LeadingDir("///"));
    EXPECT_EQ("/", mcp3d::LeadingDir("/"));
    EXPECT_EQ("/", mcp3d::LeadingDir("/foo"));
    EXPECT_EQ(".", mcp3d::LeadingDir("."));
    EXPECT_EQ(".", mcp3d::LeadingDir("./"));
    EXPECT_EQ("foo", mcp3d::LeadingDir("foo/"));
    EXPECT_EQ("foo", mcp3d::LeadingDir("foo///bar"));
    EXPECT_EQ("foo ", mcp3d::LeadingDir("foo ///bar"));
}

TEST(Utility, SplitBaseName)
{
    EXPECT_EQ(make_pair(string("/foo"), string("bar")),
              mcp3d::SplitBaseName("/foo/bar/"));
    EXPECT_EQ(make_pair(string("/foo"), string("bar")),
              mcp3d::SplitBaseName("/foo/bar//"));
    EXPECT_EQ(make_pair(string("."), string("bar")),
              mcp3d::SplitBaseName("./bar/"));
    EXPECT_EQ(make_pair(string(""), string(".")), mcp3d::SplitBaseName("./"));
    EXPECT_EQ(make_pair(string(""), string(".")), mcp3d::SplitBaseName("."));
    EXPECT_EQ(make_pair(string(""), string("/")), mcp3d::SplitBaseName("///"));
}

TEST(Utility, JoinPath)
{
    EXPECT_EQ("", mcp3d::JoinPath("", "", ""));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath({"/", "foo", "bar"}));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath({"/", "/foo/", "/bar"}));
    EXPECT_EQ("/foo/bar/", mcp3d::JoinPath({"/", "/foo/", "/bar/"}));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath({"/", "/", "/foo/", "/bar"}));
    EXPECT_EQ("./foo/bar", mcp3d::JoinPath({"./", "/foo/", "/bar"}));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath({"", "/foo/", "/bar"}));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath({"", "/", "/foo/", "/bar"}));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath("/", "foo", "bar"));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath("/", "/foo/", "/bar"));
    EXPECT_EQ("/foo/bar/", mcp3d::JoinPath("/", "/foo/", "/bar/"));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath("/", "/", "/foo/", "/bar"));
    EXPECT_EQ("./foo/bar", mcp3d::JoinPath("./", "/foo/", "/bar"));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath("", "/foo/", "/bar"));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath("", "/", "/foo/", "/bar"));
    EXPECT_EQ("/foo/bar", mcp3d::JoinPath("", "/", "/foo/", "/bar", ""));
}

TEST(Utility, PathComponents)
{
    string input("/foo/bar/");
    EXPECT_EQ(vector<string>({"/", "foo", "bar"}), mcp3d::PathComponents(input));
    input = "foo/bar/";
    EXPECT_EQ(vector<string>({"foo", "bar"}), mcp3d::PathComponents(input));
    input = "./foo///bar//";
    EXPECT_EQ(vector<string>({".", "foo", "bar"}), mcp3d::PathComponents(input));
    input = ".///";
    EXPECT_EQ(vector<string>({"."}), mcp3d::PathComponents(input));
    input = "///";
    EXPECT_EQ(vector<string>({"/"}), mcp3d::PathComponents(input));
    input = "";
    EXPECT_EQ(vector<string>({""}), mcp3d::PathComponents(input));
}

TEST(Utility, DiscardLeadingComponent)
{
    string path("/foo/bar"), component("/foo");
    EXPECT_EQ("bar", mcp3d::DiscardLeadingComponent(path, component));
    path= "/foo/bar/";
    component = "/foo";
    EXPECT_EQ("bar", mcp3d::DiscardLeadingComponent(path, component));
    path= "foo/bar/";
    component = "/foo";
    EXPECT_EQ(path, mcp3d::DiscardLeadingComponent(path, component));
    path= "///";
    component = "/";
    EXPECT_EQ("", mcp3d::DiscardLeadingComponent(path, component));
    path= "///";
    component = "foo";
    EXPECT_EQ(path, mcp3d::DiscardLeadingComponent(path, component));
    path= "///foo";
    component = "/";
    EXPECT_EQ("foo", mcp3d::DiscardLeadingComponent(path, component));
}

TEST(Utility, NearestCommonDir)
{
    vector<string> paths({"foo", "/foo", "/"});
    EXPECT_EQ(make_pair(false, string()), mcp3d::NearestCommonDir(paths));
    paths = {"/bar", "/foo", "/"};
    EXPECT_EQ(make_pair(true, string("/")), mcp3d::NearestCommonDir(paths));
    paths = {"///foo", "/foo", "/foo/"};
    EXPECT_EQ(make_pair(true, string("/foo")), mcp3d::NearestCommonDir(paths));
    paths = {"///foo", "", "/foo/"};
    EXPECT_EQ(make_pair(false, string()), mcp3d::NearestCommonDir(paths));
    paths = {"//", "///", "////"};
    EXPECT_EQ(make_pair(true, string("/")), mcp3d::NearestCommonDir(paths));
    paths = {"/foo/bar/", "/foo/bar", "/foo/bar/foobar"};
    EXPECT_EQ(make_pair(true, string("/foo/bar")), mcp3d::NearestCommonDir(paths));
    paths = {"/"};
    EXPECT_EQ(make_pair(true, string("/")), mcp3d::NearestCommonDir(paths));
    paths = {"/foo/"};
    EXPECT_EQ(make_pair(true, string("/foo")), mcp3d::NearestCommonDir(paths));
    paths = {"/foo/bar"};
    EXPECT_EQ(make_pair(true, string("/foo/bar")), mcp3d::NearestCommonDir(paths));
    paths = {""};
    EXPECT_EQ(make_pair(false, string("")), mcp3d::NearestCommonDir(paths));
}

TEST(Utility, StringLower)
{
    EXPECT_EQ("@1a b_c", mcp3d::StringLower(string("@1A B_c")));
    string s = "@1A B_c";
    mcp3d::StringLower(s);
    EXPECT_EQ("@1a b_c", s);
    EXPECT_EQ("@1a b_c", mcp3d::StringLower("@1A B_c"));
}

TEST(Utility, SplitString)
{
    string s1 = " abc def g ";
    EXPECT_TRUE(vector<string>({"abc", "def", "g"}) == mcp3d::SplitString(s1));
    EXPECT_TRUE(vector<string>({"abc", "def", "g"}) == mcp3d::SplitString(s1, " "));
    EXPECT_TRUE(vector<string>({" abc def g "}) == mcp3d::SplitString(s1, ","));
    EXPECT_TRUE(vector<string>({" abc def g "}) == mcp3d::SplitString(s1, "apple"));
    string s2 = ",abc,,def,g,";
    EXPECT_TRUE(vector<string>({"", "abc", "", "def", "g"}) == mcp3d::SplitString(s2, ","));
    EXPECT_TRUE(vector<string>({",abc,,def,g,"}) == mcp3d::SplitString(s2, " "));
    EXPECT_TRUE(vector<string>({",abc,,def,g,"}) == mcp3d::SplitString(s2, "apple"));
    string s3 = "apple is eating apple";
    EXPECT_TRUE(vector<string>({" is eating "}) == mcp3d::SplitString(s3, "apple"));
    string s4 = "appleappleapple";
    EXPECT_TRUE(vector<string>({}) == mcp3d::SplitString(s4, "apple"));
    string s5 = " appleappleapple ";
    EXPECT_TRUE(vector<string>({" ", " "}) == mcp3d::SplitString(s5, "apple"));
    string s6 = "I'm eating apple at home";
    EXPECT_TRUE(vector<string>({"I'm eating ", " at home"}) == mcp3d::SplitString(s6, "apple"));
    EXPECT_TRUE(vector<string>({"I'm eating apple at home"}) == mcp3d::SplitString(s6, "orange"));
}

TEST(Utility, PadNumStr)
{
    // invalid input number
    EXPECT_THROW(mcp3d::PadNumStr(100, 2), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::PadNumStr(0, 0), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::PadNumStr(-1, 2), mcp3d::MCPAssertionError);
    EXPECT_EQ("0", mcp3d::PadNumStr(0, 1));
    EXPECT_EQ("00", mcp3d::PadNumStr(0, 2));
    EXPECT_EQ("01", mcp3d::PadNumStr(1, 2));
    EXPECT_EQ("123", mcp3d::PadNumStr(123, 3));
    EXPECT_EQ("0123", mcp3d::PadNumStr(123, 4));
}

TEST(Utility, JoinVector)
{
    vector<int> v1 {1, 2, 3};
    EXPECT_EQ("1 2 3", mcp3d::JoinVector<int>(v1));
    EXPECT_EQ("[1, 2, 3]", mcp3d::JoinVector<int>(v1, ", ", true));
    vector<int> v2 {1};
    EXPECT_EQ("1", mcp3d::JoinVector<int>(v2));
    EXPECT_EQ("[1]", mcp3d::JoinVector<int>(v2, ", ", true));

    vector<string> s1 {"a", "bcd", "", " "};
    EXPECT_EQ("a,bcd,, ", mcp3d::JoinVector<string>(s1, ",", false));
    EXPECT_EQ("axyzbcdxyzxyz ", mcp3d::JoinVector<string>(s1, "xyz", false));
    vector<string> s2 {"abc"};
    EXPECT_EQ("abc", mcp3d::JoinVector<string>(s2, ","));
}

/// should not depend on mcp3d's filesystem functions (which are not yet tested)
TEST(Utility, SysCmdResultAndListDir)
{
    string cmd_result = mcp3d::SysCmdResult("hostname");
    EXPECT_EQ(mcp3d::HostName(), cmd_result);

    string dir = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "sys_cmd"});
    boost::filesystem::remove_all(dir);
    EXPECT_FALSE(boost::filesystem::exists(dir));
    EXPECT_TRUE(boost::filesystem::create_directories(dir));

    string subdir0 = mcp3d::JoinPath({dir, "subdir0"});
    string subdir1 = mcp3d::JoinPath({dir, "subdir1"});
    string subdir2 = mcp3d::JoinPath({dir, "subdir2"});
    EXPECT_TRUE(boost::filesystem::create_directories(subdir0));
    EXPECT_TRUE(boost::filesystem::create_directories(subdir1));
    EXPECT_TRUE(boost::filesystem::create_directories(subdir2));
    string touch_file = dir + "/a.txt";
    string touch_cmd = "touch " + touch_file;
    cmd_result = mcp3d::SysCmdResult(touch_cmd);
    EXPECT_TRUE(boost::filesystem::is_regular_file(touch_file));
    string count_dir_cmd = "ls -f " + dir + " | wc -l ";
    // should count subdir0-2, touch.txt, . and ..
    cmd_result = mcp3d::SysCmdResult(count_dir_cmd);
    EXPECT_EQ("6", cmd_result);

    vector<string> list_dir_result = mcp3d::ListDir(dir, true, false, false);
    EXPECT_EQ(vector<string>({"a.txt", "subdir0", "subdir1", "subdir2"}), list_dir_result);
    list_dir_result = mcp3d::ListDir(dir, true, true, false);
    EXPECT_EQ(vector<string>({"subdir0", "subdir1", "subdir2", "a.txt"}), list_dir_result);
    list_dir_result = mcp3d::ListDir(dir, false, true, false);
    EXPECT_EQ(vector<string>({"subdir0", "subdir1", "subdir2", "a.txt"}), list_dir_result);
    list_dir_result = mcp3d::ListDir(dir, true, true, true);
    EXPECT_EQ(vector<string>({"a.txt", "subdir2", "subdir1", "subdir0"}), list_dir_result);
    boost::filesystem::remove_all(dir);
    EXPECT_FALSE(boost::filesystem::exists(dir));
}

TEST(Utility, ResolveSymLink)
{
    string test_symlink_dir = mcp3d::JoinPath({mcp3d::test_data_dir(),
                                               "utility", "resolve_symlink"});
    if (boost::filesystem::exists(boost::filesystem::path(test_symlink_dir)))
        boost::filesystem::remove_all(boost::filesystem::path(test_symlink_dir));
    EXPECT_FALSE(boost::filesystem::exists(boost::filesystem::path(test_symlink_dir)));
    EXPECT_TRUE(boost::filesystem::create_directories(test_symlink_dir));

    // file symlink, at end
    string file_path0 = mcp3d::JoinPath({test_symlink_dir, "file0"});
    string sym_link0 = mcp3d::JoinPath({test_symlink_dir, "link0"});
    string touch_cmd = mcp3d::JoinVector<string>({"touch", file_path0}, " ", false);
    mcp3d::SysCmdResult(touch_cmd);
    string link_cmd = mcp3d::JoinVector<string>({"ln", "-s", file_path0, sym_link0}, " ", false);
    mcp3d::SysCmdResult(link_cmd);
    EXPECT_TRUE(boost::filesystem::is_symlink(boost::filesystem::path(sym_link0)));
    EXPECT_EQ(boost::filesystem::read_symlink(boost::filesystem::path(sym_link0)).string(), file_path0);
    EXPECT_EQ(file_path0, mcp3d::ResolveSymlink(sym_link0).string());

    // one directory symlink in the middle
    // symdir0 -> dir0    resolve symdir0/file1
    string d0 = mcp3d::JoinPath({test_symlink_dir, "dir0"});
    EXPECT_TRUE(boost::filesystem::create_directories(d0));
    string sym_d0 = mcp3d::JoinPath({test_symlink_dir, "sym_dir0"});
    link_cmd = mcp3d::JoinVector<string>({"ln", "-s", d0, sym_d0}, " ", false);
    mcp3d::SysCmdResult(link_cmd);
    EXPECT_TRUE(boost::filesystem::is_symlink(boost::filesystem::path(sym_d0)));
    EXPECT_EQ(boost::filesystem::read_symlink(boost::filesystem::path(sym_d0)).string(), d0);
    string file_path1 = mcp3d::JoinPath({d0, "file1"});
    string sym_link1 = mcp3d::JoinPath({sym_d0, "file1"});
    EXPECT_EQ(file_path1, mcp3d::ResolveSymlink(sym_link1).string());

    // two level directory symlinks
    // directories:
    // test_symlink_dir/subdir0
    // test_symlink_dir/subdir1/subdir2
    // regular file test_symlink_dir/subdir0/file_subdir0
    // link: test_symlink_dir/subdir0/link_subdir1 -> test_symlink_dir/subdir1
    // link: test_symlink_dir/subdir1/subdir2/link_subdir0 -> test_symlink_dir/subdir0
    // resolve test_symlink_dir/subdir0/link_subdir1/subdir2/link_subdir0/file_subdir0 = test_symlink_dir/subdir0/file_subdir0
    // make directories
    string subdir0 = mcp3d::JoinPath({test_symlink_dir, "subdir0"});
    EXPECT_TRUE(boost::filesystem::create_directories(subdir0));
    string subdir1 = mcp3d::JoinPath({test_symlink_dir, "subdir1"});
    string subdir2 = mcp3d::JoinPath({subdir1, "subdir2"});
    EXPECT_TRUE(boost::filesystem::create_directories(subdir2));
    // make file
    string file_subdir0 = mcp3d::JoinPath({subdir0, "file_subdir0"});
    touch_cmd = mcp3d::JoinVector<string>({"touch", file_subdir0}, " ", false);
    mcp3d::SysCmdResult(touch_cmd);
    EXPECT_TRUE(boost::filesystem::is_regular_file(boost::filesystem::path(file_subdir0)));
    // make links
    string link_subdir1 = mcp3d::JoinPath({subdir0, "link_subdir1"});
    boost::filesystem::create_directory_symlink(subdir1, link_subdir1);
    EXPECT_TRUE(boost::filesystem::is_symlink(link_subdir1));
    string link_subdir0 = mcp3d::JoinPath({subdir2, "link_subdir0"});
    boost::filesystem::create_directory_symlink(subdir0, link_subdir0);
    EXPECT_TRUE(boost::filesystem::is_symlink(link_subdir0));
    string start = mcp3d::JoinPath({subdir0, "link_subdir1", "subdir2", "link_subdir0", "file_subdir0"}),
           resolute = mcp3d::JoinPath({subdir0, "file_subdir0"});
    // boost::filesystem::read_symlink follow link once
    EXPECT_EQ(subdir1, boost::filesystem::read_symlink(link_subdir1));
    EXPECT_TRUE(boost::filesystem::is_symlink(link_subdir1));
    EXPECT_EQ(subdir0, boost::filesystem::read_symlink(link_subdir0));
    EXPECT_EQ(resolute, mcp3d::ResolveSymlink(start).string());
    boost::filesystem::remove_all(test_symlink_dir);
    EXPECT_FALSE(boost::filesystem::exists(test_symlink_dir));
}

TEST(Utility, IsFile)
{
    string file_path = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility",
                                        "isfile_test"}),
           symlink_path = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility",
                                           "isfile_test_symlink0"}),
           dir_path = mcp3d::JoinPath(mcp3d::test_data_dir(), "utility",
                                      "isfile_test_dir");
    // return true on file path
    string touch_command = mcp3d::JoinVector<string>({"touch", file_path}, " ", false);
    mcp3d::SysCmdResult(touch_command);
    EXPECT_TRUE(boost::filesystem::is_regular_file(file_path));
    EXPECT_TRUE(boost::filesystem::is_regular_file(file_path));
    EXPECT_FALSE(boost::filesystem::is_symlink(file_path));
    // return true on symlink to file path
    boost::filesystem::create_symlink(file_path, symlink_path);
    EXPECT_TRUE(boost::filesystem::is_regular_file(symlink_path));
    EXPECT_TRUE(boost::filesystem::is_symlink(symlink_path));
    EXPECT_TRUE(mcp3d::IsFile(symlink_path));
    // return false on directory path
    EXPECT_TRUE(boost::filesystem::create_directories(dir_path));
    EXPECT_TRUE(boost::filesystem::is_directory(dir_path));
    EXPECT_FALSE(mcp3d::IsFile(dir_path));
    // clean up
    boost::filesystem::remove(file_path);
    boost::filesystem::remove(symlink_path);
    boost::filesystem::remove_all(dir_path);
    EXPECT_FALSE(mcp3d::IsFile(file_path));
    EXPECT_FALSE(boost::filesystem::exists(file_path));
    EXPECT_FALSE(mcp3d::IsFile(symlink_path));
    EXPECT_FALSE(boost::filesystem::exists(symlink_path));
    EXPECT_FALSE(mcp3d::IsFile(dir_path));
    EXPECT_FALSE(boost::filesystem::exists(dir_path));
}

TEST(Utility, MakeDirectories)
{
    string dir0 = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "make_directories"});
    string dir1 = mcp3d::JoinPath({dir0, "more", "levels", "of", "directories"});
    EXPECT_FALSE(boost::filesystem::is_directory(dir0));
    mcp3d::MakeDirectories(dir1);
    EXPECT_TRUE(boost::filesystem::is_directory(dir1));
    boost::filesystem::remove_all(dir0);
    EXPECT_FALSE(boost::filesystem::is_directory(dir0));
}

void PopulateRemovePathDir(const string& dir_path)
{
    mcp3d::MakeDirectories(dir_path);
    EXPECT_TRUE(boost::filesystem::is_directory(dir_path));
    string file_path1(mcp3d::JoinPath(dir_path, "f1")),
           file_path2(mcp3d::JoinPath(dir_path, "f2"));
    ofstream f1(file_path1);
    f1.close();
    ofstream f2(file_path2);
    f2.close();
    EXPECT_TRUE(boost::filesystem::is_regular_file(file_path2));
    EXPECT_TRUE(boost::filesystem::is_regular_file(file_path2));
}

TEST(Utility, RemovePath)
{
    // remove regular file and directories
    string dir = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "remove_path"});
    PopulateRemovePathDir(dir);
    EXPECT_TRUE(mcp3d::RemovePath(mcp3d::JoinPath({dir, "f1"})));
    EXPECT_FALSE(boost::filesystem::is_regular_file(mcp3d::JoinPath({dir, "f1"})));
    EXPECT_TRUE(mcp3d::RemovePath(dir));
    EXPECT_FALSE(boost::filesystem::is_regular_file(mcp3d::JoinPath({dir, "f2"})));
    EXPECT_FALSE(boost::filesystem::is_directory(dir));

    // giving non existing path returns true
    EXPECT_TRUE(mcp3d::RemovePath(mcp3d::JoinPath({mcp3d::test_data_dir(), "does_not_exist"})));

    // multi-threading remove regular file and directories
    dir = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "multi_thread_remove_path"});
    PopulateRemovePathDir(dir);
    int n_threads = 4;
    #pragma omp parallel num_threads(n_threads)
    {
        cout << "remove path from thread " << omp_get_thread_num() << endl;
        EXPECT_TRUE(mcp3d::RemovePath(dir));
    }

    // remove symlink with/without its target
    dir = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "remove_link_and_target"});
    PopulateRemovePathDir(dir);
    // make chain of symlink link0 -> link1 -> dir
    string link0 = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "link0"}),
           link1 = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "link1"});
    boost::filesystem::create_directory_symlink(dir, link1);
    boost::filesystem::create_symlink(link1, link0);
    EXPECT_TRUE(boost::filesystem::is_symlink(link0));
    EXPECT_EQ(dir, boost::filesystem::read_symlink(link1).string());
    EXPECT_TRUE(boost::filesystem::is_directory(link0));
    EXPECT_TRUE(boost::filesystem::is_symlink(link1));
    EXPECT_EQ(link1, boost::filesystem::read_symlink(link0).string());
    EXPECT_TRUE(boost::filesystem::is_directory(link1));
    EXPECT_TRUE(boost::filesystem::is_directory(dir));
    // remove link0 without remove its target
    EXPECT_TRUE(mcp3d::RemovePath(link0, false));
    EXPECT_FALSE(boost::filesystem::exists(link0));
    EXPECT_TRUE(boost::filesystem::is_symlink(link1));
    EXPECT_TRUE(boost::filesystem::is_directory(dir));
    // remove link0 with its target
    boost::filesystem::create_symlink(link1, link0);
    EXPECT_TRUE(boost::filesystem::is_symlink(link0));
    EXPECT_EQ(link1, boost::filesystem::read_symlink(link0).string());
    EXPECT_TRUE(boost::filesystem::is_directory(link0));
    EXPECT_TRUE(mcp3d::RemovePath(link0));
    EXPECT_FALSE(boost::filesystem::exists(link0));
    EXPECT_FALSE(boost::filesystem::exists(link1));
    EXPECT_FALSE(boost::filesystem::exists(dir));

    // multi-threading remove symlink with/without its target
    dir = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "multi_thread_remove_link_and_target"});
    PopulateRemovePathDir(dir);
    // make chain of symlink link0 -> link1 -> dir
    boost::filesystem::create_directory_symlink(dir, link1);
    boost::filesystem::create_symlink(link1, link0);
    EXPECT_TRUE(boost::filesystem::is_symlink(link0));
    EXPECT_EQ(dir, boost::filesystem::read_symlink(link1).string());
    EXPECT_TRUE(boost::filesystem::is_directory(link0));
    EXPECT_TRUE(boost::filesystem::is_symlink(link1));
    EXPECT_EQ(link1, boost::filesystem::read_symlink(link0).string());
    EXPECT_TRUE(boost::filesystem::is_directory(link1));
    EXPECT_TRUE(boost::filesystem::is_directory(dir));
    #pragma omp parallel num_threads(n_threads)
    {
        cout << "remove link and targets from thread " << omp_get_thread_num() << endl;
        EXPECT_TRUE(mcp3d::RemovePath(link0));
    }
    EXPECT_FALSE(boost::filesystem::exists(link0));
    EXPECT_FALSE(boost::filesystem::exists(link1));
    EXPECT_FALSE(boost::filesystem::exists(dir));
}

TEST(Utility, EmptyDir)
{
    string dir = mcp3d::JoinPath(mcp3d::test_data_dir(), "empty_dir");

}

TEST(Utility, FilesInDir)
{
    string dir = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility", "files_in_dir"});
    EXPECT_TRUE(boost::filesystem::is_directory(dir));
    vector<string> include({".tif", ".ims"});
    vector<string> exclude({"~"});

    bool full_path = false, sort_files = true;
    vector<string> all_files = mcp3d::FilesInDir(dir, full_path, sort_files);
    EXPECT_EQ(vector<string>({"a.ims", "b.tif", "c~.tif", "d.jpg"}), all_files);
    vector<string> incl_files = mcp3d::FilesInDir(dir, full_path, sort_files, include);
    EXPECT_EQ(vector<string>({"a.ims", "b.tif", "c~.tif"}), incl_files);
    vector<string> incl_2files = mcp3d::FilesInDir(dir, full_path, sort_files, include, 2);
    EXPECT_EQ(vector<string>({"a.ims", "b.tif"}), incl_2files);
    vector<string> incl_excl_files = mcp3d::FilesInDir(dir, full_path, sort_files, include, -1, exclude);
    EXPECT_EQ(vector<string>({"a.ims", "b.tif"}), incl_excl_files);

    full_path = true;
    vector<string> incl_excl_full_files = mcp3d::FilesInDir(dir, full_path, sort_files, include, 1, exclude);
    EXPECT_EQ(vector<string>({mcp3d::JoinPath(dir, "a.ims")}), incl_excl_full_files);
}

TEST(Utility, DirsInDir)
{
    string dir = mcp3d::JoinPath({mcp3d::test_data_dir(), "utility",
                                  "dirs_in_dir"});
    vector<string> dir_names {"pyr_ratio_1", "pyr", "apyrratio", "random", "pyr_ratio_exclude"},
                   include {"pyr", "ratio"}, exclude {"exclude"},
                   dir_names_incl {"pyr_ratio_1", "apyrratio", "pyr_ratio_exclude"},
                   dir_names_incl_excl {"pyr_ratio_1", "apyrratio"};
    for (const auto& dir_name: dir_names)
        EXPECT_TRUE(mcp3d::MakeDirectories(mcp3d::JoinPath(dir, dir_name)));

    vector<string> result, expected;
    for (const auto& dir_name: dir_names)
        expected.push_back(mcp3d::JoinPath(dir, dir_name));
    sort(expected.begin(), expected.end());
    EXPECT_EQ(expected, mcp3d::DirsInDir(dir, true, true));
    expected.pop_back();
    expected.pop_back();
    EXPECT_EQ(expected, mcp3d::DirsInDir(dir, true, true, {}, 3, {}));
    expected.clear();

    for (const auto& dir_name: dir_names)
        expected.push_back(dir_name);
    sort(expected.begin(), expected.end());
    EXPECT_EQ(expected, mcp3d::DirsInDir(dir, false, true));
    expected.clear();

    for (const auto& dir_name: dir_names_incl)
        expected.push_back(mcp3d::JoinPath(dir, dir_name));
    sort(expected.begin(), expected.end());
    EXPECT_EQ(expected, mcp3d::DirsInDir(dir, true, true, include));
    expected.clear();

    for (const auto& dir_name: dir_names_incl)
        expected.push_back(dir_name);
    sort(expected.begin(), expected.end());
    EXPECT_EQ(expected, mcp3d::DirsInDir(dir, false, true, include));
    expected.clear();

    for (const auto& dir_name: dir_names_incl_excl)
        expected.push_back(mcp3d::JoinPath(dir, dir_name));
    sort(expected.begin(), expected.end());
    EXPECT_EQ(expected, mcp3d::DirsInDir(dir, true, true, include, -1, exclude));
    expected.clear();

    for (const auto& dir_name: dir_names_incl_excl)
        expected.push_back(dir_name);
    sort(expected.begin(), expected.end());
    EXPECT_EQ(expected, mcp3d::DirsInDir(dir, false, true, include, -1, exclude));
    expected.clear();

    EXPECT_TRUE(mcp3d::RemovePath(dir));
}

TEST(Utility, IntPow)
{
    EXPECT_THROW(mcp3d::IntPow(-1, 0), mcp3d::MCPDomainError);
    EXPECT_THROW(mcp3d::IntPow(0, -1), mcp3d::MCPDomainError);
    EXPECT_EQ(0, mcp3d::IntPow(0, 2));
    EXPECT_EQ(1, mcp3d::IntPow(0, 0));
    EXPECT_EQ(1, mcp3d::IntPow(2, 0));
    EXPECT_EQ(2, mcp3d::IntPow(2, 1));
    EXPECT_EQ(64, mcp3d::IntPow(2, 6));
}

TEST(Utility, DFSDirectoryWalker)
{
    /*  walker test dir
     *  ---- 1
     *       ---- 11
     *             ---- 111.txt 112.txt 113.txt
     *       ---- 12
     *       ---- 13
     *  ---- 2
     *       ---- 21
     *              ---- 211
     *                     ---- 2111
     *                             ---- 2111.txt 2112.txt 2113.txt
     * ---- 3
     *      ---- 31 32 33
     * ---- 1.txt  2.txt 3.txt
     */
    string walker_test_dir = mcp3d::JoinPath(mcp3d::test_data_dir(),
                                             "utility", "dfs_directory_walker");
    EXPECT_TRUE(mcp3d::IsDir(walker_test_dir));
    bool sort_alpha = true;
    mcp3d::DFSDirectoryWalker walker(walker_test_dir, sort_alpha);
    EXPECT_TRUE(walker.TerminalFiles().empty());
    EXPECT_EQ(walker_test_dir, walker.TerminalDirPath());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir/1/11, files 111.txt 112.txt 113.txt, n_dirs = 0
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"111.txt", "112.txt", "113.txt"}, string(1, '\n')), walker.TerminalFiles());
    EXPECT_EQ(mcp3d::JoinPath({walker_test_dir, "1/11"}), walker.TerminalDirPath());
    EXPECT_EQ(0, walker.TerminalDirCounts());
    EXPECT_EQ(3, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir/1, files 11 12 13, n_dirs = 3
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"11", "12", "13"}, string(1, '\n')), walker.TerminalFiles());
    EXPECT_EQ(mcp3d::JoinPath({walker_test_dir, "1"}), walker.TerminalDirPath());
    EXPECT_EQ(3, walker.TerminalDirCounts());
    EXPECT_EQ(3, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir/2/21/211/2111, files 2111.txt 2112.txt 2113.txt, n_dirs = 0
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"2111.txt", "2112.txt", "2113.txt"}, string(1, '\n')), walker.TerminalFiles());
    EXPECT_EQ(mcp3d::JoinPath({walker_test_dir, "2/21/211/2111"}), walker.TerminalDirPath());
    EXPECT_EQ(0, walker.TerminalDirCounts());
    EXPECT_EQ(3, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir/2/21/211, files 2111, n_dirs = 1
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"2111"}), walker.TerminalFiles());
    EXPECT_EQ(mcp3d::JoinPath({walker_test_dir, "2/21/211"}), walker.TerminalDirPath());
    EXPECT_EQ(1, walker.TerminalDirCounts());
    EXPECT_EQ(1, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir/2/21, files 211, n_dirs = 1
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"211"}), walker.TerminalFiles());
    EXPECT_EQ(mcp3d::JoinPath({walker_test_dir, "2/21"}), walker.TerminalDirPath());
    EXPECT_EQ(1, walker.TerminalDirCounts());
    EXPECT_EQ(1, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir/2, files 21, n_dirs = 1
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"21"}), walker.TerminalFiles());
    EXPECT_EQ(mcp3d::JoinPath({walker_test_dir, "2"}), walker.TerminalDirPath());
    EXPECT_EQ(1, walker.TerminalDirCounts());
    EXPECT_EQ(1, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir/3, files 31 32 33, n_dirs = 3
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"31", "32", "33"}, string(1, '\n')), walker.TerminalFiles());
    EXPECT_EQ(mcp3d::JoinPath({walker_test_dir, "3"}), walker.TerminalDirPath());
    EXPECT_EQ(3, walker.TerminalDirCounts());
    EXPECT_EQ(3, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker_test_dir, files 1 2 3 1.txt 2.txt 3.txt, n_dirs = 3
    walker.ExpandNextTerminal();
    EXPECT_FALSE(walker.TerminalFiles().empty());
    EXPECT_EQ(mcp3d::JoinVector<string>({"1", "2", "3", "1.txt", "2.txt", "3.txt"}, string(1, '\n')), walker.TerminalFiles());
    EXPECT_EQ(walker_test_dir, walker.TerminalDirPath());
    EXPECT_EQ(3, walker.TerminalDirCounts());
    EXPECT_EQ(6, walker.TerminalEntryCounts());
    EXPECT_FALSE(walker.Exausted());

    // -> walker exausted
    walker.ExpandNextTerminal();
    EXPECT_TRUE(walker.TerminalFiles().empty());
    EXPECT_EQ("", walker.TerminalDirPath());
    EXPECT_TRUE(walker.Exausted());

    // no error should be thrown by attempts to expand exhausted tree
    EXPECT_NO_THROW(walker.ExpandNextTerminal());
    EXPECT_TRUE(walker.TerminalFiles().empty());
    EXPECT_EQ("", walker.TerminalDirPath());
    EXPECT_TRUE(walker.Exausted());

    // empty directory case
    string empty_walker_test_dir = mcp3d::JoinPath(mcp3d::test_data_dir(),
                                             "utility", "empty_dfs_directory_walker");
    mcp3d::MakeDirectories(empty_walker_test_dir);
    mcp3d::DFSDirectoryWalker empty_walker(empty_walker_test_dir, sort_alpha);
    // immediate exhaustion
    empty_walker.ExpandNextTerminal();
    EXPECT_TRUE(walker.TerminalFiles().empty());
    EXPECT_EQ("", empty_walker.TerminalDirPath());
    EXPECT_TRUE(walker.Exausted());

    // file path case
    string file_walker_test_path = mcp3d::JoinPath(mcp3d::test_data_dir(),
                                                   "utility", "file_dfs_directory_walker");
    mcp3d::SysCmdResult(mcp3d::JoinVector<string>({"touch", file_walker_test_path}, " ", false));
    EXPECT_TRUE(mcp3d::IsFile(file_walker_test_path));

    mcp3d::DFSDirectoryWalker file_walker(file_walker_test_path, sort_alpha);
    // immediate exhaustion
    empty_walker.ExpandNextTerminal();
    EXPECT_TRUE(walker.TerminalFiles().empty());
    EXPECT_EQ("", empty_walker.TerminalDirPath());
    EXPECT_TRUE(walker.Exausted());
    boost::filesystem::remove(file_walker_test_path);
}

TEST(Utility, AsyncRemovePath)
{
    // create random directory structures
    string robo_remove_dir = mcp3d::JoinPath(mcp3d::test_data_dir(),
                                             "utility", "async_remove");
    long seed;
    INIT_TIMER(0)
    cout << "creating 5 random directory structures" << endl;
    for (int i = 0; i < 5; ++i)
    {
        seed = chrono::system_clock::now().time_since_epoch().count();

        // non blocking deletion
        mcp3d::MakeDirectories(robo_remove_dir);
        EXPECT_TRUE(boost::filesystem::is_directory(robo_remove_dir));
        TIC(0)
        int n_items = mcp3d::test::CreateDirectoryTree(robo_remove_dir, seed);
        cout << "directory " << robo_remove_dir << " contains " << n_items << " items" << endl;
        TOC(0)
        REPORT_TIME_TO_COMPLETION("time to create directory structure", 0)
        TIC(0)
        mcp3d::AsyncRemovePath(robo_remove_dir);
        TOC(0)
        REPORT_TIME_TO_COMPLETION("time to remove directory structure with async io", 0)
        EXPECT_TRUE(!boost::filesystem::is_directory(robo_remove_dir));

        // blocking deletion
        mcp3d::MakeDirectories(robo_remove_dir);
        EXPECT_TRUE(boost::filesystem::is_directory(robo_remove_dir));
        TIC(0)
        n_items = mcp3d::test::CreateDirectoryTree(robo_remove_dir, seed);
        cout << "directory " << robo_remove_dir << " contains " << n_items << " items" << endl;
        TOC(0)
        REPORT_TIME_TO_COMPLETION("time to create directory structure", 0)
        TIC(0)
        mcp3d::SysCmdResult(mcp3d::JoinVector<string>({"rm", "-r", robo_remove_dir}, " "));
        TOC(0)
        REPORT_TIME_TO_COMPLETION("time to remove directory structure with blocking io", 0)
        EXPECT_TRUE(!boost::filesystem::is_directory(robo_remove_dir));
        cout << endl;
    }
}

TEST(Utility, AllGreater)
{
    vector<int> v0({10, 11, 23}), v1({10, 10, 22});
    EXPECT_FALSE(mcp3d::AllGreater(v0, v1));
    vector<int> v2({10, 11, 23}), v3({9, 10, 22});
    EXPECT_TRUE(mcp3d::AllGreater(v2, v3));
}

TEST(Utility, To5D)
{
    vector<int> dims0({10, 10}), dims1({3, 3, 3, 3, 3}),
                empty_dims, long_dims({1, 1, 1, 10, 10, 1});
    EXPECT_EQ(vector<int>({1, 1, 1, 10, 10}), mcp3d::To5D(dims0));
    EXPECT_EQ(dims1, mcp3d::To5D(dims1));
    EXPECT_THROW(mcp3d::To5D(empty_dims), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::To5D(long_dims), mcp3d::MCPAssertionError);

}

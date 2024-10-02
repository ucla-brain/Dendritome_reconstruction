//
// Created by muyezhu on 5/24/19.
//
#include "common/mcp3d_utility.hpp"
#include "image/mcp3d_image_utils.hpp"
#include "app2_parameters.hpp"

using namespace std;

string mcp3d::App2Parameters::MetaString()
{
    stringstream meta_stream;
    meta_stream << "# foreground_percent = " << foreground_percent_ << endl;
    meta_stream << "# foreground_percent_by_plane = " << foreground_percent_by_plane_ << endl;
    meta_stream << "# threshold_value = " << threshold_value_ << endl;
    meta_stream << "# plane_threshold values = " << plane_threshold_values_string() << endl;
    meta_stream << "# length_thresh = " << length_thresh_ << endl;
    meta_stream << "# SR_ratio = " << sr_ratio_ << endl;
    meta_stream << "# gsdt = " << gsdt_ << endl;
    meta_stream << "# allow_gap = " << allow_gap_ << endl;
    meta_stream << "# cnn_type = " << cnn_type_ << endl;
    meta_stream << "# radius_from_2d = " << radius_from_2d_  << endl;
    meta_stream << "# swc_resample = " << swc_resample_  << endl;
    meta_stream << "# cube_256 = " << cube_256_  << endl;
    meta_stream << "# marker_file_path = " << (marker_file_path_.empty() ? "none" : marker_file_path_)  << endl;
    return meta_stream.str();
}

void mcp3d::App2CommandLineArgs::PrintUsage()
{
    cout<< "Usage : ./app2 <volume_root_dir> <channel> [-somda-dir <soma_dir>] [-inmarker <marker_file>] [-outswc <swc_file>] "
            "[-resolution-level <int>] [-image-offsets <int> [<int>] [<int>]] [-image-extents <int> [<int>] [<int>]]"
            "[-thresh-value <int>] [-fg-percent <double>] [-by-plane] [-gsdt] [-cnn-type <int>] [-length-thresh <double>] [-allow-gap]" << endl;
    cout << endl;
    cout << "-soma-dir           [-sd] soma directory. if not given it's assumed images under volume_root_dir has filled soma" << endl;
    cout << "-inmarker           [-im] input marker file, the first marker is source marker, the rest are target markers" << endl;
    cout << "-outswc             [-os] output tracing result, default is <imagename>_tracing.swc" << endl;
    cout << "-resolution-level   [-rl] resolution level to perform tracing at. default is 0, ie original resolution" << endl;
    cout << "-image-offsets      [-io] offsets of subvolume, in zyx order. each axis has default offset value 0 if not provided" << endl;
    cout << "-image-extents      [-ie] extents of subvolume, in zyx order. each axis has extent from offset to axis maximum range if not provided" << endl;
    cout << "-thresh-value       [-tv] threshold value used in GSDT and tree construction when no target marker" << endl;
    cout << "-fg-percent         [-fp] percent of voxels to be considered foreground. if positive, overrides -thresh-val" << endl;
    cout << "-by-plane           [-bp] perform foreground percent threshold by each z plane. ignored if -thresh-val is given and -fg-percent is not" << endl;
    cout << "-length-thresh      [-lt] default 1.0, the length threshold value for hierarchy pruning(hp)"<<endl;
    cout << "-sr-ratio           [-sr] default 1/3, signal/redundancy ratio threshold" << endl;
    cout << "-gsdt               [-gs] perform GSDT for original image" << endl;
    cout << "-cnn-type           [-ct] default 2. connection type 1 for 6 neighbors, 2 for 18 neighbors, 3 for 26 neighbors" << endl;
    cout << "-allow-gap          [-ag] accept one background point between foreground during tree construction when only one marker" << endl;
    cout << endl;
    cout << "Hierarchical pruning method by Hang Xiao"<<endl;
    cout << endl;
}

string mcp3d::App2CommandLineArgs::MetaString()
{
    stringstream meta_stream;
    meta_stream << "# image root dir = " << image_root_dir_ << endl;
    meta_stream << "# channel = " << channel_number_ << endl;
    meta_stream << "# resolution level = " << resolution_level_ << endl;
    meta_stream << "# offsets (zyx) = " << mcp3d::JoinVector(image_offsets(), ", ", true) << endl;
    meta_stream << "# extents (zyx) = " << mcp3d::JoinVector(image_extents(), ", ", true) << endl;
    meta_stream << app2_parameters_.MetaString();
    return meta_stream.str();
}

bool mcp3d::ParseApp2Args(int argc, char * argv[], mcp3d::App2CommandLineArgs& args)
{
    if (argc < 3)
    {
        mcp3d::App2CommandLineArgs::PrintUsage();
        return false;
    }
    try
    {
        // global volume and channel selection
        args.set_image_root_dir(argv[1]);
        args.set_channel_number(ParseInt(argv[2]));
        // if the switch is given, parameter(s) corresponding to the switch is expected
        for (int i = 3; i < argc; ++i)
        {
            // subvolume selection arguments
            if (strcmp(argv[i],"-resolution-level") == 0 || strcmp(argv[i],"-rl") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -resolution-level")
                args.set_resolution_level(ParseInt(argv[i + 1]));
                ++i;
            }
            else if (strcmp(argv[i],"-image-offsets") == 0 || strcmp(argv[i],"-io") == 0)
            {
                vector<int> offsets;
                for (int j = 0; j < 3; ++ j)
                {
                    if (i + 1 >= argc || argv[i + 1][0] == '-')
                        MCP3D_INVALID_ARGUMENT("missing parameter(s) for -image-offsets")
                    int offset = ParseInt(argv[i + 1]);
                    offsets.push_back(offset);
                    ++i;
                }
                args.set_image_offsets(offsets);
            }
            else if (strcmp(argv[i],"-image-extents") == 0 || strcmp(argv[i],"-ie") == 0)
            {
                vector<int> extents;
                for (int j = 0; j < 3; ++ j)
                {
                    if (i + 1 >= argc || argv[i + 1][0] == '-')
                        MCP3D_INVALID_ARGUMENT("missing parameter for -image-extents")
                    int extent = ParseInt(argv[i + 1]);
                    extents.push_back(extent);
                    ++i;
                }
                args.set_image_extents(extents);
            }
            else if (strcmp(argv[i], "-outswc") == 0 || strcmp(argv[i], "-os") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -outswc")
                args.set_swc_path(argv[i + 1]);
                ++i;
            }
            else if (strcmp(argv[i], "-inmarker") == 0 || strcmp(argv[i], "-im") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -outswc")
                args.app2_parameters().set_marker_file_path(argv[i + 1]);
                ++i;
            }
            else if (strcmp(argv[i], "-thresh-value") == 0 || strcmp(argv[i], "-tv") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -thresh-value")
                args.app2_parameters().set_threshold_value(ParseInt(argv[i + 1]));
                ++i;
            }
            else if (strcmp(argv[i], "-fg-percent") == 0 || strcmp(argv[i], "-fp") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -fg-percent")
                args.app2_parameters().set_foreground_percent(ParseDouble(argv[i + 1]));
                ++i;
            }
            else if (strcmp(argv[i], "-by-plane") == 0 || strcmp(argv[i], "-bp") == 0)
                args.app2_parameters().set_foreground_percent_by_plane(true);
            else if (strcmp(argv[i], "-length-thresh") == 0 || strcmp(argv[i], "-lt") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -length-thresh")
                args.app2_parameters().set_length_thresh(ParseDouble(argv[i + 1]));
                ++i;
            }
            else if (strcmp(argv[i], "-sr-ratio") == 0 || strcmp(argv[i], "-sr") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -sr-ratio")
                args.app2_parameters().set_sr_ratio(ParseDouble(argv[i + 1]));
                ++i;
            }
            else if (strcmp(argv[i], "-cnn-type") == 0 || strcmp(argv[i], "-ct") == 0)
            {
                if (i + 1 >= argc || argv[i + 1][0] == '-')
                    MCP3D_INVALID_ARGUMENT("missing parameter for -cnn-type")
                args.app2_parameters().set_cnn_type(ParseInt(argv[i + 1]));
                ++i;
            }
            else if (strcmp(argv[i], "-gsdt") == 0 || strcmp(argv[i], "-gs") == 0)
                args.app2_parameters().set_gsdt(true);
            else if (strcmp(argv[i], "-allow-gap") == 0 || strcmp(argv[i], "-ag") == 0)
                args.app2_parameters().set_allow_gap(true);
            else
                cout << "unknown option " << argv[i] << endl;
        }
        // give default swc path if not given
        if (args.swc_path().empty())
        {
            string z_start = to_string(args.image_offsets()[0]),
                   y_start = to_string(args.image_offsets()[1]),
                   x_start = to_string(args.image_offsets()[2]);
            string z_end = to_string(args.image_offsets()[0] + args.image_extents()[0]),
                   y_end = to_string(args.image_offsets()[1] + args.image_extents()[1]),
                   x_end = to_string(args.image_offsets()[2] + args.image_extents()[2]);
            args.set_swc_path(mcp3d::JoinPath(args.image_root_dir(),
                                              "tracing_z" + z_start + "_" + z_end +
                                              "_y" + y_start + "_" + y_end +
                                              "_x" + x_start + "_" + x_end) + ".swc");
        }
        // if output directory for swc file does not exist, create it
        mcp3d::MakeDirectories(mcp3d::ParentDir(args.swc_path()));
        // if neither background threshold nor foreground percent given, set foreground percent to 0.01
        if (args.app2_parameters().threshold_value() < 0 && args.app2_parameters().foreground_percent() < 0)
            args.app2_parameters().set_foreground_percent(0.01);
        return true;
    }
    catch (const exception& e)
    {
        cout << e.what() << endl;
        MCP3D_MESSAGE("invalid command line arguments. neuron tracing not performed")
        App2CommandLineArgs::PrintUsage();
        return false;
    }
}

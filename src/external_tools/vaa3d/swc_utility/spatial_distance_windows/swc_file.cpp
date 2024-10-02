//
// Created by muyezhu on 3/3/19.
//
#include <cstdint>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include "swc_file.hpp"

using namespace std;


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

vector<string> SplitString(const string &in_string, const string &delim)
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

string StripLeading(const std::string &input)
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

string StripTrailing(const std::string &input)
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

string Strip(const std::string &input)
{
    string trim_leading = StripLeading(input);
    return StripTrailing(trim_leading);
}

string ReadMetaLines(const string& swc_path)
{
    ifstream input(swc_path, ifstream::in);
    if (!input.good())
    {
        cout << "can not open " << swc_path << endl;
        return string();
    }
    string line = string(), trimmed_line, meta_lines;
    while(getline(input, line))
    {
        trimmed_line = StripLeading(line);
        if (trimmed_line.empty())
            continue;
        if (trimmed_line[0] == '#')
        {
            meta_lines.append(trimmed_line);
            meta_lines.append("\n");
        }
        else
            break;
    }
    return meta_lines;
}

NeuronTree readSWC_file_to_tree(const string& file_name)
{
    NeuronTree nt;
    nt.file_name = file_name;
    ifstream input_file(file_name, ifstream::in);
    if (!input_file.good())
    {
        cout << "can not open " << file_name << endl;
        return nt;
    }

    int count = 0;
    vector<NeuronSWC> listNeuron;
    unordered_map<int64_t , int64_t>  hashNeuron;
    listNeuron.clear();
    hashNeuron.clear();

    string line, trimmed_line;
    while (getline(input_file, line))
    {
        trimmed_line = Strip(line);
        if (trimmed_line.empty())
            continue;
        if (trimmed_line[0]=='\0' || trimmed_line[0] == '#')
            continue;

        count++;
        NeuronSWC S;

        vector<string> tokens = SplitString(trimmed_line, " ");
        if (tokens.empty())
            continue;
        try
        {
            for (size_t i=0; i<tokens.size(); i++)
            {
                if (i==0) S.n = stoi(tokens[i]);
                else if (i==1) S.type = stoi(tokens[i]);
                else if (i==2) S.x = stod(tokens[i]);
                else if (i==3) S.y = stod(tokens[i]);
                else if (i==4) S.z = stod(tokens[i]);
                else if (i==5) S.r = stod(tokens[i]);
                else if (i==6) S.pn = stoi(tokens[i]);
                    //the ESWC extension, by PHC, 20120217
                else if (i==7) S.seg_id = stoi(tokens[i]);
                else if (i==8) S.level = stoi(tokens[i]);
                    //change ESWC format to adapt to flexible feature number, by WYN, 20150602
                else
                    S.fea_val.push_back(stof(tokens[i]));
            }
        }
        catch (const exception& e)
        {
            cout << "bad swc line: " + trimmed_line << endl;
        }

        listNeuron.push_back(S);
        hashNeuron.insert(make_pair(S.n, listNeuron.size() - 1));

    }

    if (listNeuron.size()<1)
        return nt;
    //now update other NeuronTree members
    nt.n = 1; //only one neuron if read from a file
    nt.listNeuron = listNeuron;
    nt.hashNeuron = hashNeuron;
    return nt;
}

bool writeSWC_file(const string& filename, const NeuronTree& nt, const string& infostring)
{
    ofstream output_file(filename, ofstream::out);
    if (!output_file.good())
    {
        cout << "could not create " + filename << endl;
        return false;
    }

    if (!infostring.empty())
        output_file << infostring;

    NeuronSWC * p_pt=0;
    char buffer[1024];
    for (int i=0;i<nt.listNeuron.size(); i++)
    {
        p_pt = (NeuronSWC *)(&(nt.listNeuron.at(i)));
        sprintf(buffer, "%ld %d %5.3f %5.3f %5.3f %5.3f %ld\n",
                p_pt->n, p_pt->type, p_pt->x, p_pt->y, p_pt->z, p_pt->r, p_pt->pn);
        output_file << buffer;
    }
    output_file.close();
}

vector<MyMarker*> readSWC_file(string swc_file)
{
    vector<MyMarker*> swc;

    ifstream ifs(swc_file.c_str());

    if(ifs.fail())
    {
        cout << "open swc file : " << swc_file << " error" << endl;
        return swc;
    }

    map<int, MyMarker*> marker_map;
    map<MyMarker*, int> parid_map;
    while(ifs.good())
    {
        if(ifs.peek() == '#'){ifs.ignore(1000,'\n'); continue;}
        MyMarker *  marker = new MyMarker;
        int my_id = -1 ; ifs >> my_id;
        if(my_id == -1) break;
        if(marker_map.find(my_id) != marker_map.end())
        {
            cerr<<"Duplicate Node. This is a graph file. Please read is as a graph."<<endl; //return vector<MyMarker*>();
        }
        marker_map[my_id] = marker;

        ifs>> marker->type;
        ifs>> marker->x;
        ifs>> marker->y;
        ifs>> marker->z;
        ifs>> marker->radius;
        int par_id = -1; ifs >> par_id;

        parid_map[marker] = par_id;
        swc.push_back(marker);
    }
    ifs.close();
    vector<MyMarker*>::iterator it = swc.begin();
    while(it != swc.end())
    {
        MyMarker * marker = *it;
        marker->parent = marker_map[parid_map[marker]];
        it++;
    }
    return swc;
}

bool readSWC_file(string swc_file, vector<MyMarker> & outmarkers)
{
    ifstream ifs(swc_file.c_str());

    if(ifs.fail())
    {
        cout<<"open swc file : "<< swc_file <<" error"<<endl;
        return false;
    }

    while(ifs.good())
    {
        if(ifs.peek() == '#'){ifs.ignore(1000,'\n'); continue;}
        MyMarker  marker;
        int my_id = -1 ; ifs >> my_id;
        if(my_id == -1) break;

        ifs>> marker.type;
        ifs>> marker.x;
        ifs>> marker.y;
        ifs>> marker.z;
        ifs>> marker.radius;
        int par_id = -1; ifs >> par_id;
        //cout<<"("<<marker.x<<","<<marker.y<<","<<marker.z<<")"<<endl;

        outmarkers.push_back(marker);
    }
    ifs.close();
    return true;
}

bool saveSWC_file(const string &swc_path, vector<MyMarker *> &outmarkers, const string& meta_lines)
{
    list<string> swc_info_string;
    if (!meta_lines.empty())
        swc_info_string.push_back(meta_lines);
    return saveSWC_file(swc_path, outmarkers, swc_info_string);
}

bool saveSWC_file(const string& swc_path, vector<MyMarker*> & outmarkers, list<string> & info_strings)
{
    if(swc_path.find_last_of(".dot") == swc_path.size() - 1)
        return saveDot_file(swc_path, outmarkers);

    cout <<"marker num = " << outmarkers.size() << ", save swc file to " << swc_path << endl;
    map<MyMarker*, int> ind;
    ofstream ofs(swc_path.c_str());

    if(ofs.fail())
    {
        cout<<"open swc file error"<<endl;
        return false;
    }

    list<string>::iterator it;
    for (const auto& info_string: info_strings)
        ofs << info_string << endl;

    ofs<<"##n,type,x,y,z,radius,parent"<<endl;
    for(int i = 0; i < outmarkers.size(); i++) ind[outmarkers[i]] = i+1;

    for(int i = 0; i < outmarkers.size(); i++)
    {
        MyMarker * marker = outmarkers[i];
        int parent_id;
        if(marker->parent == 0) parent_id = -1;
        else parent_id = ind[marker->parent];
        ofs << i + 1 << " " << marker->type << " "
            << marker->x << " " <<marker->y << " " << marker->z
            << " " <<marker->radius << " " << parent_id<<endl;
    }
    ofs.close();
    return true;
}

bool saveSWC_file(const string& swc_path, vector<NeuronSWC*> & outmarkers, list<string> & meta_lines)
{
    cout<<"marker num = "<<outmarkers.size()<<", save swc file to "<<swc_path<<endl;
    ofstream ofs(swc_path.c_str());

    if(ofs.fail())
    {
        cout<<"open swc file error"<<endl;
        return false;
    }
    ofs<<"#name "<<swc_path<<endl;
    ofs<<"#comment "<<endl;

    list<string>::iterator it;
    for (it=meta_lines.begin();it!=meta_lines.end(); it++)
        ofs<< *it <<endl;

    ofs<<"##n,type,x,y,z,radius,parent"<<endl;

    for(int64_t i = 0; i < outmarkers.size(); i++)
    {
        NeuronSWC * marker = outmarkers[i];
        int64_t parent_id = (marker->parent < 0) ? -1 : marker->parent;
        ofs<<i<<" "<<marker->type<<" "<<marker->x<<" "<<marker->y<<" "<<marker->z<<" "<<marker->radius<<" "<<parent_id<<endl;
    }
    ofs.close();
    return true;
}

bool saveSWC_file(const string& swc_path, const NeuronTree& neuron_tree, const string& meta_lines)
{
    ofstream out_file(swc_path, fstream::out);
    if (!out_file.good())
    {
        cout << "can not open file " << swc_path << " for write" << endl;
        return false;
    }

    if (!meta_lines.empty())
        out_file << meta_lines << endl;

    for (size_t i = 0; i < neuron_tree.listNeuron.size(); i++)
    {
        out_file << neuron_tree.listNeuron.at(i).n << " "
                 << neuron_tree.listNeuron.at(i).type << " "
                 << neuron_tree.listNeuron.at(i).x << " "
                 << neuron_tree.listNeuron.at(i).y << " "
                 << neuron_tree.listNeuron.at(i).z << " "
                 << neuron_tree.listNeuron.at(i).r << " "
                 << neuron_tree.listNeuron.at(i).pn << endl;
    }
    out_file.close();
    return true;
}

bool saveDot_file(const string& swc_path, vector<MyMarker*> & outmarkers)
{
    cout<<"marker num = "<<outmarkers.size()<<", save swc file to "<<swc_path<<endl;
    map<MyMarker*, int> ind;
    ofstream ofs(swc_path.c_str());

    if(ofs.fail())
    {
        cout<<"open swc file error"<<endl;
        return false;
    }
    ofs<<"digraph \""<<swc_path<<"\" {"<<endl;
    ofs<<"\trankdir = BT;"<<endl;

    for(int i = 0; i < outmarkers.size(); i++) ind[outmarkers[i]] = i+1;
    for(int i = 0; i < outmarkers.size(); i++)
    {
        MyMarker * marker = outmarkers[i];
        if(marker->parent)
        {
            int parent_id = ind[marker->parent];
            MyMarker * parent = marker->parent;
            ofs<<"\t"<<i+1<<" -> "<<parent_id<<";"<<endl;
        }
    }
    ofs<<"}"<<endl;
    ofs.close();
    return true;
}


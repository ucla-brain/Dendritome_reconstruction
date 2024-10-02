//
// Created by muyezhu on 5/29/19.
//
#include <map>
#include "common/mcp3d_macros.hpp"
#include "vaa3d/swc_utility/swc_file.hpp"
#include "resample_swc.hpp"

using namespace std;

#define DISTP(a,b) sqrt(((a)->x-(b)->x)*((a)->x-(b)->x)+((a)->y-(b)->y)*((a)->y-(b)->y)+((a)->z-(b)->z)*((a)->z-(b)->z))

struct Point
{
    double x,y,z,r;
    int64_t type;
    Point* p;
    int64_t childNum;
    int64_t level,seg_id;
    vector<float> fea_val;
};

typedef vector<Point*> Segment;
typedef vector<Point*> Tree;

void resample_path(Segment * seg, double step)
{
    char c;
    Segment seg_r;
    double path_length = 0;
    Point* start = seg->at(0);
    Point* seg_par = seg->back()->p;
    int64_t iter_old = 0;
    seg_r.push_back(start);
    while (iter_old < seg->size() && start && start->p)
    {
        path_length += DISTP(start,start->p);
        if (path_length<=seg_r.size()*step)
        {
            start = start->p;
            iter_old++;
        }
        else//a new point should be created
        {
            path_length -= DISTP(start,start->p);
            Point* pt = new Point;
            double rate = (seg_r.size()*step-path_length)/(DISTP(start,start->p));
            pt->x = start->x + rate*(start->p->x-start->x);
            pt->y = start->y + rate*(start->p->y-start->y);
            pt->z = start->z + rate*(start->p->z-start->z);
            pt->r = start->r*(1-rate) + start->p->r*rate;//intepolate the radius
            pt->p = start->p;

            if (rate<0.5)
            {
                pt->type = start->type;
                pt->seg_id = start->seg_id;
                pt->level = start->level;
                pt->fea_val = start->fea_val;
            }
            else
            {
                pt->type = start->p->type;
                pt->seg_id = start->p->seg_id;
                pt->level = start->p->level;
                pt->fea_val = start->p->fea_val;

            }
            seg_r.back()->p = pt;
            seg_r.push_back(pt);
            path_length += DISTP(start,pt);
            start = pt;
        }
    }
    seg_r.back()->p = seg_par;
    for (int64_t i=0;i<seg->size();i++)
        if (!seg->at(i)) {delete seg->at(i); seg->at(i) = NULL;}
    *seg = seg_r;
};

NeuronTree ResampleNeuronTree(NeuronTree input, double step)
{
    NeuronTree result;
    int64_t siz = input.listNeuron.size();
    Tree tree;
    for (int64_t i=0;i<siz;i++)
    {
        NeuronSWC s = input.listNeuron[i];
        Point* pt = new Point;
        pt->x = s.x;
        pt->y = s.y;
        pt->z = s.z;
        pt->r = s.r;
        pt ->type = s.type;
        pt->seg_id = s.seg_id;
        pt->level = s.level;
        pt->fea_val = s.fea_val;
        pt->p = NULL;
        pt->childNum = 0;
        tree.push_back(pt);
    }
    for (int64_t i=0;i<siz;i++)
    {
        if (input.listNeuron[i].pn<0)
            continue;
        int64_t pid = input.hashNeuron.at(input.listNeuron[i].pn);
        tree[i]->p = tree[pid];
        tree[pid]->childNum++;
    }
//	printf("tree constructed.\n");
    vector<Segment*> seg_list;
    for (int64_t i=0;i<siz;i++)
    {
        if (tree[i]->childNum!=1)//tip or branch point
        {
            Segment* seg = new Segment;
            Point* cur = tree[i];
            do
            {
                seg->push_back(cur);
                cur = cur->p;
            }
            while(cur && cur->childNum==1);
            seg_list.push_back(seg);
        }
    }
//	printf("segment list constructed.\n");
    for (int64_t i=0;i<seg_list.size();i++)
    {
        resample_path(seg_list[i], step);
    }

//	printf("resample done.\n");
    tree.clear();
    map<Point*, int64_t> index_map;
    for (int64_t i=0;i<seg_list.size();i++)
        for (int64_t j=0;j<seg_list[i]->size();j++)
        {
            tree.push_back(seg_list[i]->at(j));
            index_map.insert(pair<Point*, int64_t>(seg_list[i]->at(j), tree.size()-1));
        }
    for (int64_t i=0;i<tree.size();i++)
    {
        NeuronSWC S;
        Point* p = tree[i];
        S.n = i+1;
        if (p->p==NULL) S.pn = -1;
        else
            S.pn = index_map[p->p]+1;
        if (p->p==p) printf("There is loop in the tree!\n");
        S.x = p->x;
        S.y = p->y;
        S.z = p->z;
        S.r = p->r;
        S.type = p->type;
        S.seg_id = p->seg_id;
        S.level = p->level;
        S.fea_val = p->fea_val;
        result.listNeuron.push_back(S);
    }
    for (int64_t i=0;i<tree.size();i++)
    {
        if (tree[i]) {delete tree[i]; tree[i]=NULL;}
    }
    for (int64_t j=0;j<seg_list.size();j++)
        if (seg_list[j]) {delete seg_list[j]; seg_list[j] = NULL;}
    for (int64_t i=0;i<result.listNeuron.size();i++)
        result.hashNeuron.insert(make_pair(result.listNeuron[i].n, i));
    return result;
}

void mcp3d::ResampleSwc(const string& input_swc_path, double step, const string& output_swc_path)
{
    if (step <= 0)
    {
        MCP3D_MESSAGE("expecting step to be positive but received negative value. setting to 10.0")
        step = 10.0;
    }
    string output_path = output_swc_path.empty()?
                         input_swc_path.substr(0, input_swc_path.size() - 4) + "_resample.swc" : output_swc_path;
    NeuronTree nt = readSWC_file_to_tree(input_swc_path);
    NeuronTree result = ResampleNeuronTree(nt, step);
    string meta_lines = mcp3d::ReadMetaLines(input_swc_path);
    meta_lines.append("# resampled from ");
    meta_lines.append(input_swc_path);
    meta_lines.append("\n# resample step = " + to_string(step));
    saveSWC_file(output_path, result, meta_lines);
    cout << "resampled swc saved to " << output_path << endl;
}


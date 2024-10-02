/*
 * Copyright (c)2006-2010  Hanchuan Peng (Janelia Farm, Howard Hughes Medical Institute).  
 * All rights reserved.
 */


/************
                                            ********* LICENSE NOTICE ************

This folder contains all source codes for the V3D project, which is subject to the following conditions if you want to use it. 

You will ***have to agree*** the following terms, *before* downloading/using/running/editing/changing any portion of codes in this package.

1. This package is free for non-profit research, but needs a special license for any commercial purpose. Please contact Hanchuan Peng for details.

2. You agree to appropriately cite this work in your related studies and publications.

Peng, H., Ruan, Z., Long, F., Simpson, J.H., and Myers, E.W. (2010) “V3D enables real-time 3D visualization and quantitative analysis of large-scale biological image data sets,” Nature Biotechnology, Vol. 28, No. 4, pp. 348-353, DOI: 10.1038/nbt.1612. ( http://penglab.janelia.org/papersall/docpdf/2010_NBT_V3D.pdf )

Peng, H, Ruan, Z., Atasoy, D., and Sternson, S. (2010) “Automatic reconstruction of 3D neuron structures using a graph-augmented deformable model,” Bioinformatics, Vol. 26, pp. i38-i46, 2010. ( http://penglab.janelia.org/papersall/docpdf/2010_Bioinfo_GD_ISMB2010.pdf )

3. This software is provided by the copyright holders (Hanchuan Peng), Howard Hughes Medical Institute, Janelia Farm Research Campus, and contributors "as is" and any express or implied warranties, including, but not limited to, any implied warranties of merchantability, non-infringement, or fitness for a particular purpose are disclaimed. In no event shall the copyright owner, Howard Hughes Medical Institute, Janelia Farm Research Campus, or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; reasonable royalties; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.

4. Neither the name of the Howard Hughes Medical Institute, Janelia Farm Research Campus, nor Hanchuan Peng, may be used to endorse or promote products derived from this software without specific prior written permission.

*************/




//by Hanchuan Peng
//2009, Jan
//last update: 090206
//last update 091127. write a new joining program
//last update: 2011-05-12 fix a small potential bug

#ifndef MCP3D_VAA3D_SD_V_NEURON_SWC_HPP
#define MCP3D_VAA3D_SD_V_NEURON_SWC_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <map>


struct V_NeuronSWC_coord    //for sort
{
	union {
	double data[3];
	struct { double x,y,z; };
	};
	bool equal(const V_NeuronSWC_coord & other) const { return x==other.x && y==other.y && z==other.z;}
	bool equal(double x1, double y1, double z1) const { return x==x1 && y==y1 && z==z1;}
	bool nonequal(const V_NeuronSWC_coord &other) const {return !equal(other);}
	bool nonequal(double x1, double y1, double z1) const {return !equal(x1, y1, z1);}
	void set(double x1, double y1, double z1) {x=x1; y=y1; z=z1;}
};

inline bool operator == (const V_NeuronSWC_coord a, const V_NeuronSWC_coord b)
{
	return (a.x==b.x && a.y==b.y && a.z==b.z);
}

inline bool operator < (const V_NeuronSWC_coord a, const V_NeuronSWC_coord b)
{
	return ((a.x<b.x) ||
			(a.x==b.x && a.y<b.y) ||
			(a.x==b.x && a.y==b.y && a.z<b.z)
			);
}

inline double distL2square(const V_NeuronSWC_coord & a, const V_NeuronSWC_coord & b)
{
	return ((a.x-b.x)*(a.x-b.x) + 
			(a.y-b.y)*(a.y-b.y) + 
			(a.z-b.z)*(a.z-b.z) );
}

struct V_NeuronSWC_unit
{
	union {
		double data[10];
		struct {
			double n, type, x, y, z, r, parent,
			nchild,
			seg_id, nodeinseg_id;
		};
	};
    V_NeuronSWC_unit() {for (int64_t i=0;i<int64_t(sizeof(data)/sizeof(double));i++) data[i]=0; r=0.5;}
	operator V_NeuronSWC_coord() {V_NeuronSWC_coord c; c.x=x; c.y=y; c.z=z; return c;}
	V_NeuronSWC_coord get_coord() {V_NeuronSWC_coord c; c.x=x; c.y=y; c.z=z; return c;}
	void set(double x1, double y1, double z1, double r1, double p1, double t1) {x=x1; y=y1; z=z1; r=r1;parent=p1;type=t1;}
	void set(double x1, double y1, double z1, double r1, double p1) {x=x1; y=y1;z=z1;r=r1;parent=p1;}
	void set(double x1, double y1, double z1, double r1) {x=x1; y=y1;z=z1;r=r1;}
	void set(double x1, double y1, double z1) {x=x1; y=y1;z=z1;}
};

inline double distL2square(const V_NeuronSWC_unit & a, const V_NeuronSWC_unit & b)
{
	return ((a.x-b.x)*(a.x-b.x) + 
			(a.y-b.y)*(a.y-b.y) + 
			(a.z-b.z)*(a.z-b.z) );
}

struct Node_Link  // for graph form of swc
{
	int64_t i;
	std::vector<int64_t> in_link;
    std::vector<int64_t> out_link;
	int64_t nlink;
};
typedef std::map<int64_t, Node_Link> Link_Map;
typedef std::map<int64_t, Node_Link>::iterator Link_Map_Iter;

struct V_SWCNodes
{
    std::vector<int64_t> nid_array; //the array of node id (node id is the 1st column value)
    std::vector<V_NeuronSWC_coord> ncoord_array; //following the same order of nid_array
    std::map<int64_t, int64_t> nid_ipos_lut; //a look up table to return the index of an nid, i.e. nid_ipos_lut[nid_array.at(i)] equals i;
};

struct V_NeuronSWC
{
    std::vector<V_NeuronSWC_unit> row;
	bool b_linegraph;
    std::string name;
    std::string comment;
    std::string file;
	unsigned char color_uc[4];
	bool b_jointed;
    bool to_be_deleted;   // @ADDED by Alessandro on 2015-05-08. Needed to support late delete of multiple neuron segments.
    bool to_be_broken;
	bool on; //Added by Y. Wang on 2016-05-25. For the segment-wise display of a SWC.

	bool check_data_consistency() {/* to verify if unique node id have unique coord, and if parent are in the nid, except -1*/ return true;}

	V_NeuronSWC(std::string new_name="unset", bool is_linegraph=false)
	{
		name=new_name; b_linegraph=is_linegraph;  *(int*)color_uc=0; b_jointed=false;
        to_be_deleted = false;
        to_be_broken = false;
		on = true;
	}

	void printInfo();

	int64_t nrows() {return row.size();}

	V_SWCNodes unique_nodes_info(); //determine the unique nodes

    std::vector<int64_t> unique_nid(); //unique node id (id is the first column value in SWC)
	int64_t n_unique_nid(); //number of unique node ids
    std::vector<V_NeuronSWC_coord> unique_ncoord(); //unique node coordinates (coordinates are the 3rd to 5th column)
	int64_t n_unique_ncoord(); //number of unique node coords

	int64_t maxnoden() //091029 change maxnoden from >=-1 to >=0 for base_n in set_simple_path...
	{
                int64_t maxn=0;	for (int64_t i=0;i<(int64_t)row.size();i++) if (row.at(i).n > maxn) maxn = row.at(i).n;		return maxn;
	}
	int64_t getIndexofParent(int64_t j)
	{
		int64_t res=-1; int64_t parent = row.at(j).parent;
                for (int64_t i=0;i<(int64_t)row.size();i++) if (row.at(i).n==parent)	{res=i; break;}
		return res;
	}
	std::vector<int64_t> getIndexofParent_nodeid(int64_t nid) //return the array of of node "nid"'s parents' nid
	{
		std::vector<int64_t> res;
                for (int64_t i=0;i<(int64_t)row.size();i++)
		{
			if (row.at(i).n==nid)
			{
				int64_t curparent = row.at(i).parent;
				bool b_exist=false;
                                for (int64_t j=0;j<(int64_t)res.size();j++)
					if (res.at(j)==curparent) {	b_exist=true; break;}
				if (!b_exist)
					res.push_back(curparent);
			}
		}
		return res;
	}

	void append(V_NeuronSWC_unit & new_node) {row.push_back(new_node);}
	void clear() {row.clear();}
	std::vector<V_NeuronSWC> decompose();
	bool reverse();

	bool isLineGraph() {return b_linegraph;} //just return the "claimed" property is a line graph
	//check if a 3D location is contained in the swc
	int64_t getFirstIndexof3DPos(double x,double y,double z) //return -1 is no included, othwise return the first detected index
	{
		int64_t res=-1;
                for (int64_t i=0;i<(int64_t)row.size();i++) if (row.at(i).data[2]==x && row.at(i).data[3]==y && row.at(i).data[4]==z)	{res=i; break;}
		return res;
	}
	int64_t getFirstIndexof3DPos(const V_NeuronSWC_unit & subject_node) {return getFirstIndexof3DPos(subject_node.data[2], subject_node.data[3], subject_node.data[4]);}
	int64_t getFirstIndexof3DPos(const V_NeuronSWC_unit * subject_node) {return getFirstIndexof3DPos(subject_node->data[2], subject_node->data[3], subject_node->data[4]);}

	std::vector<int64_t> getAllIndexof3DPos(double x,double y,double z, int64_t noninclude_ind) //return all indexes except the one indicated as noninclude_ind
	{
		std::vector<int64_t> res;
                for (int64_t i=0;i<(int64_t)row.size();i++) if (row.at(i).data[2]==x && row.at(i).data[3]==y && row.at(i).data[4]==z)	{ if (i!=noninclude_ind) res.push_back(i); }
		return res;
	}
	std::vector<int64_t> getAllIndexof3DPos(const V_NeuronSWC_unit & subject_node, int64_t noninclude_ind) {return getAllIndexof3DPos(subject_node.data[2], subject_node.data[3], subject_node.data[4], noninclude_ind);}
	std::vector<int64_t> getAllIndexof3DPos(const V_NeuronSWC_unit * subject_node, int64_t noninclude_ind) {return getAllIndexof3DPos(subject_node->data[2], subject_node->data[3], subject_node->data[4], noninclude_ind);}
};

struct V_NeuronSWC_list
{
	std::vector<V_NeuronSWC> seg; //since each seg could be a complete neuron or multiple paths, thus I call it "seg", but not "path"
	int64_t last_seg_num; //?? for what purpose? seems only used once in v3d_core.cpp. Questioned by Hanchuan, 20100210
	std::string name;
	std::string comment;
	std::string file;
	unsigned char color_uc[4];
	bool b_traced;

	V_NeuronSWC_list() {last_seg_num=-1; *(int*)color_uc=0; b_traced=true;}

	int64_t nsegs() {return seg.size();}
        int64_t nrows() {int64_t n=0; for (int64_t i=0;i<(int64_t)seg.size();i++) n+=seg.at(i).nrows(); return n;}
	int64_t maxnoden()
	{
                int64_t maxn=0;	for (int64_t i=0;i<(int64_t)seg.size();i++) if (seg.at(i).maxnoden() > maxn) maxn = seg.at(i).maxnoden();	return maxn;
	}
	bool isJointed() {return nsegs()==1 && seg.at(0).b_jointed;}

	void append(V_NeuronSWC & new_seg) {seg.push_back(new_seg); last_seg_num=seg.size();}
        void append(std::vector<V_NeuronSWC> & new_segs) {for (int k=0; k<(int)new_segs.size(); k++) seg.push_back(new_segs.at(k)); last_seg_num=seg.size();}
	void clear() {last_seg_num=seg.size(); seg.clear();}
	void merge();
	void decompose();
	bool reverse();
	bool split(int64_t seg_id, int64_t nodeinseg_id);
    bool deleteSeg(int64_t seg_id);

    // @ADDED by Alessandro on 2015-05-08. Needed to support late delete of multiple neuron segments.
    void                                            // no value returned
        deleteMultiSeg(                             // by default, deletes neuron segments having 'to_be_deleted' field set to 'true'
            std::vector <int64_t> *seg_ids = 0);    // if provided, deletes the corresponding neuron segments.
};

bool verifyIsLineGraph(const V_NeuronSWC & in_swc); //this will use graph algorithm to verify if really a line graph as claimed

///////////////////////////////
//// NOTE: merge_xxx() and join_xxx() operations are different! merge() means just putting together, but join_xxx() means more, also including removing redundant nodes.
//////////////////////////////

Link_Map get_link_map(const V_NeuronSWC & in_swc);
std::vector<V_NeuronSWC> decompose_V_NeuronSWC(V_NeuronSWC & in_swc);
V_NeuronSWC join_V_NeuronSWC_vec(std::vector<V_NeuronSWC> & in_swc_vec);

bool reverse_V_NeuronSWC_inplace(V_NeuronSWC & in_swc);
bool change_type_in_seg_of_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id, int type);
bool change_radius_in_seg_of_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id, double radius);

V_NeuronSWC merge_V_NeuronSWC_list(V_NeuronSWC_list & in_swc_list);
bool delete_seg_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id); //delete a seg in the V_NeuronSWC_list
double length_seg_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id); //find the length of a seg
double getLength_V_NeuronSWC(V_NeuronSWC & subject_swc); //compute the length of a swc neuron
V_NeuronSWC join_segs_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id_array[], int n_segs); //merge several segs in V_NeuronSWC_list, indexed by seg_id_array (length is n_segs)
bool join_two_V_NeuronSWC_old(V_NeuronSWC & destination_swc, V_NeuronSWC & subject_swc); //add the subject_swc neuron to the destination_swc, and merge all overlaping nodes
bool join_two_V_NeuronSWC(V_NeuronSWC & destination_swc, V_NeuronSWC & subject_swc); //add the subject_swc neuron to the destination_swc, and merge all overlaping nodes

int64_t find_node_in_V_NeuronSWC(V_NeuronSWC & in_swc, double x, double y, double z); //find the id of a node
int64_t find_seg_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, double x, double y, double z, int64_t & nodeinseg_id ); //find the id of a seg

V_NeuronSWC_list split_V_NeuronSWC_simplepath(V_NeuronSWC & in_swc, int64_t nodeinseg_id);
V_NeuronSWC_list split_V_NeuronSWC_simplepath(V_NeuronSWC & in_swc, double x, double y, double z);

std::map<int64_t,int64_t> unique_V_NeuronSWC_nodeindex(V_NeuronSWC & in_swc);
bool simplify_V_NeuronSWC_nodeindex(V_NeuronSWC & my_swc); //map the node index of a swc neuron to the range 1~N (N is the number of unique indexes)



template <class T> //should be a struct at least included members (x,y,z)
void set_simple_path_unit (V_NeuronSWC_unit & v, int64_t base_n, std::vector<T> & mUnit, int64_t i, bool link_order, double r=1, double default_type=3)
{
	int64_t N = mUnit.size();
		v.type	= default_type;
		v.x 	= mUnit[i].x;
		v.y 	= mUnit[i].y;
		v.z 	= mUnit[i].z;
		v.r 	= r;
	if (link_order) // same as index order
	{
		v.n		= base_n +1+i;
		v.parent = (i>=N-1)? -1 : (v.n +1);
	}
	else // reversed link order
	{
		v.n		= base_n +1+i;
		v.parent = (i<=0)? -1 : (v.n -1);
	}
}
template <class T> //should be a struct at least included members (x,y,z)
void set_simple_path (V_NeuronSWC & cur_seg, int64_t base_n, std::vector<T> & mUnit, bool link_order, double r=1, double default_type=3)
{
	cur_seg.clear();
	for (int64_t i=0;i<mUnit.size();i++)
	{
		V_NeuronSWC_unit v;
        set_simple_path_unit(v, base_n, mUnit, i, link_order, r, default_type);
		cur_seg.append(v);
		//qDebug("%d ", cur_seg.nnodes());
	}
}
template <class T> //should be a struct at least included members (n, parent)
void reset_simple_path_index (int64_t base_n, std::vector<T> & mUnit)
{
	int64_t N = mUnit.size();
	for (int64_t i=0; i<mUnit.size(); i++)
	{
		if (mUnit[0].parent >=1) // same as index order
		{
			mUnit[i].n = base_n +1+i;
			mUnit[i].parent = (i>=N-1)? -1: (mUnit[i].n +1);
		}
		else                    // reversed link order
		{
			mUnit[i].n = base_n +1+i;
			mUnit[i].parent = (i<=0)? -1: (mUnit[i].n -1);
		}
	}
}
template <class T>
int64_t value_in_vector(std::vector<T> & vec, T & val) //a <0 return value indicating the value does not exist in vector
{
	if (vec.size()<=0) return -1;
	for (int64_t i=0;i<vec.size();i++)
		if (vec.at(i)==val)
			return i;
	return -1;
}
template <class T>
std::vector<T> unique_values_in_vector(std::vector<T> & vec) //return a new vector where value are unique value in vec
{
	if (vec.size()<=0) return vec;
	std::vector<T> vnew;
	vnew.push_back(vec.at(0));
	for (int64_t i=1;i<vec.size();i++)
	{
		if (value_in_vector(vnew, vec.at(i))<0) //not exist yet
			vnew.push_back(vec.at(i));
	}
	return vnew;
}


#endif //MCP3D_VAA3D_SD_V_NEURON_SWC_HPP



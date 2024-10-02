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
//last update: 090430
//last update 091127. add a new neuron_joining function

#include <cstdio>
#include <cmath>
#include <iostream>
#include "common/mcp3d_macros.hpp"
#include "v_neuronswc.hpp"

using namespace std;

//V_NeuronSWC:: /////////////////////////////////////////////////
void V_NeuronSWC::printInfo()
{
	MCP3D_MESSAGE("name=[" + name + "] comment=[" + comment + "] file=[" + file +
                  "] b_linegraph=[" + to_string(b_linegraph) +
                  "] b_jointed=[" + to_string(b_jointed) + "]")
	for (int64_t i=0;i<row.size();i++)
	{
		V_NeuronSWC_unit v = row.at(i);
		cout << "row [" << i << "] : " << v.data[0] << " " << v.data[1]
             << " " << v.data[2] << " " << v.data[3] << " " << v.data[4]
             << " " << v.data[5] << " " << v.data[6] << endl;
	}
    cout << "*********** finish print info" << endl;
}


vector<V_NeuronSWC> V_NeuronSWC::decompose()
{
	return decompose_V_NeuronSWC(*this);
}

bool V_NeuronSWC::reverse()
{
	return reverse_V_NeuronSWC_inplace(*this);
}

//V_NeuronSWC_list:: ////////////////////////////////////////////
void V_NeuronSWC_list::merge()
{
	V_NeuronSWC new_seg = merge_V_NeuronSWC_list(*this);
	clear();
	append(new_seg);
}

void V_NeuronSWC_list::decompose()
{
	vector <V_NeuronSWC> new_segs;
	new_segs.clear();
	for (int k=0; k<seg.size(); k++)
	{
		vector <V_NeuronSWC> tmp_segs = seg.at(k).decompose();
		for (int j=0; j<tmp_segs.size(); j++)
		{
			new_segs.push_back(tmp_segs.at(j));
		}
	}
	clear();
	append(new_segs);
}

bool V_NeuronSWC_list::reverse()
{
	int oldnseg = seg.size();
	bool res = false;
	//must first decompose to simple segments
	{
		seg = merge_V_NeuronSWC_list(*this).decompose();
	}
	for (int i=0; i<seg.size(); i++)
	{
		res = reverse_V_NeuronSWC_inplace(seg[i]);
		if (!res)  break;
	}
	// keep old segment number
	if (oldnseg==1)
	{
		merge();
	}
	return res;
}

bool V_NeuronSWC_list::split(int64_t seg_id, int64_t nodeinseg_id)
{
	bool res = false;
	V_NeuronSWC_list new_slist = split_V_NeuronSWC_simplepath (seg.at(seg_id), nodeinseg_id);

	if (!delete_seg_in_V_NeuronSWC_list(*this, seg_id))
	{
		printf("error in deleting seg (%d)\n", seg_id);
		res = false;
	}
	else
	{
		append(new_slist.seg);
		res = true;
	}
	return res;
}

bool V_NeuronSWC_list::deleteSeg(int64_t seg_id)
{
	return delete_seg_in_V_NeuronSWC_list(*this, seg_id);
}

// @ADDED by Alessandro on 2015-05-08. Needed to support late delete of multiple neuron segments.
void                                            // no value returned
    V_NeuronSWC_list::deleteMultiSeg(           // by default, deletes neuron segments having 'to_be_deleted' field set to 'true'
        std::vector <int64_t> *seg_ids /*= 0*/) // if provided, deletes the corresponding neuron segments.
{
    // first pass (optional) to mark segments to delete
    if(seg_ids)
        for(int k=0; k<seg_ids->size(); k++)
            if((*seg_ids)[k] >= 0 && (*seg_ids)[k]<seg.size())
                seg[(*seg_ids)[k]].to_be_deleted = true;

    // second pass to delete segments
    std::vector<V_NeuronSWC>::iterator iter = seg.begin();
    while (iter != seg.end())
        if (iter->to_be_deleted)
            iter = seg.erase(iter);
        else
           ++iter;

    last_seg_num=seg.size();
}

/////////////////////////////////////////////////////////////////////////////////

V_NeuronSWC merge_V_NeuronSWC_list(V_NeuronSWC_list & in_swc_list)
{
	V_NeuronSWC out_swc;  out_swc.clear();
	V_NeuronSWC_unit v;
	int64_t n=0, i,j,k;
	int64_t nsegs = in_swc_list.seg.size();
	for (k=0;k<nsegs;k++)
    {
        if(in_swc_list.seg.at(k).to_be_deleted)
            continue;

		if(!in_swc_list.seg.at(k).on)
            continue;

		vector <V_NeuronSWC_unit> &row = (in_swc_list.seg.at(k).row);
		if (row.size()<=0) continue;

		//first find the min index number, then all index will be automatically adjusted
		int64_t min_ind = row[0].n;
		for (j=1;j<row.size();j++)
		{
			if (row[j].n < min_ind)  min_ind = row[j].n;
			if (min_ind<0)
                cout << "Found illeagal neuron node index which is less than 0 in merge_V_NeuronSWC_list()!" << endl;
		}
		//qDebug()<<min_ind;

		// segment id & color type
		int seg_id = k;

		//now merge
		int64_t n0=n;
		for (j=0;j<row.size();j++)
		{
            v.seg_id = seg_id;
            v.nodeinseg_id = j;

			v.n = (n0+1) + row[j].n-min_ind;
			for (i=1;i<=5;i++)	v.data[i] = row[j].data[i];
			v.parent = (row[j].parent<0)? -1 : ((n0+1) + row[j].parent-min_ind); //change row[j].parent<=0 to row[j].parent<0, PHC 091123.

			//qDebug()<<row[j].n<<"->"<<v.n<<" "<<row[j].parent<<"->"<<v.parent<<" "<<n;

			out_swc.row.push_back(v);
            n++;
		}
	}
	out_swc.color_uc[0] = in_swc_list.color_uc[0];
	out_swc.color_uc[1] = in_swc_list.color_uc[1];
	out_swc.color_uc[2] = in_swc_list.color_uc[2];
	out_swc.color_uc[3] = in_swc_list.color_uc[3];

	out_swc.name = "merged";
//	for (i=0;i<out_swc.nrows();i++)
//		qDebug()<<out_swc.row.at(i).data[2]<<" "<<out_swc.row.at(i).data[3]<<" "<<out_swc.row.at(i).data[4]<<" "<<out_swc.row.at(i).data[6];

	return out_swc;
}


//090925 RZC
vector <V_NeuronSWC> decompose_V_NeuronSWC_old(V_NeuronSWC & in_swc)
{
	// map swc's index --> vector's index
	map<double,int64_t> index_map;
	index_map.clear();
	for (int64_t i=0; i<in_swc.row.size(); i++)
	{
		index_map[in_swc.row.at(i).n] = i;
		in_swc.row[i].nchild = 0;
	}

	// count parent->nchild
	for (int64_t i=0; i<in_swc.row.size(); i++)
	{
		double parent = in_swc.row.at(i).parent;
		if (parent>0) // 1-based
		{
			int64_t j = index_map[parent]; // this is fast
			in_swc.row[j].nchild ++;
		}
	}

	vector <V_NeuronSWC> out_swc_segs;
	out_swc_segs.clear();

	for (;;)
	{
		// check is all nodes processed
		int64_t n_removed = 0;
		for (int64_t i=0; i<in_swc.row.size(); i++)
		{
			V_NeuronSWC_unit & cur_node = in_swc.row[i];
			if (cur_node.nchild<0) // is removed node
				n_removed ++;
		}
		if (n_removed >= in_swc.row.size())
		{
			cout << "split_V_NeuronSWC_segs is done." << endl;
			break; //over, all nodes have been labeled to remove
		}

		// find a tip point as start point
        int64_t j = -1;
		int n_left = 0;
		for (int64_t i=0; i < in_swc.row.size(); i++)
		{
			V_NeuronSWC_unit & cur_node = in_swc.row[i];
			if (cur_node.nchild ==0) // tip point, assume link from leaf to root
			{
				j = i;
				break; //find
			}
			else if (cur_node.nchild >=1) // path point
			{
				n_left++;
			}
		}
		if (j<0) //not find
		{
			if (n_left)
			{
				cout << "split_V_NeuronSWC_segs cann't find tip point, the SWC has loop link!" << endl;
			}
			break;
		}

		// extract a simple segment
		V_NeuronSWC new_seg;
		new_seg.clear();
		//qDebug("decompose_V_NeuronSWC_segs: segment from node #%d", j);

		for (int64_t n=1; ;n++)
		{
			V_NeuronSWC_unit & cur_node = in_swc.row[j];
			cur_node.nchild --; // processed 1 time

			V_NeuronSWC_unit new_node = cur_node;
			new_node.n = n;
			new_node.parent = n+1; // same link order as original
			new_seg.row.push_back(new_node);

			if(cur_node.parent <1)    // root point ////////////////////////////
			{
				//qDebug("decompose_V_NeuronSWC_segs: segment end at root #%d", j);
				cur_node.nchild = floor(cur_node.nchild) + 0.5; //0.5 to prevent from treating as tip point
				break; //over, a simple segment
			}
			else if (cur_node.nchild >0)  // branch point (after --) ///////////
			{
				//qDebug("decompose_V_NeuronSWC_segs: segment end at branch #%d", j);
				cur_node.nchild = floor(cur_node.nchild) + 0.5; //0.5 to prevent from treating as tip point
				if (cur_node.nchild <1)  cur_node.nchild = 0;   //all children processed, set it as tip point
				break; //over, a simple segment
			}
			else    // path node ///////////////////////////////////////////////
			{
				//qDebug("decompose_V_NeuronSWC_segs: node #%d", j);
				double parent = cur_node.parent;
				j = index_map[parent]; //// this is fast
				cur_node.nchild = -1;  // label to remove path node
			}
		}

		if (new_seg.row.size()>=2)//==0? single point
		{
			new_seg.row[new_seg.row.size()-1].parent = -1; // set segment end
			char buf[10];
			new_seg.name = sprintf(buf,"%d", out_swc_segs.size()+1);
			new_seg.b_linegraph=true; //donot forget to do this //####################
			out_swc_segs.push_back(new_seg);
		}
	}

	return out_swc_segs;
}

//091212 RZC
Link_Map get_link_map(const V_NeuronSWC & in_swc)
{
	// map swc's index --> vector's index & in/out link
	Link_Map link_map;
	link_map.clear();

	// setting index of map
	for (int64_t i=0; i<in_swc.row.size(); i++)
	{
        int64_t nid = in_swc.row.at(i).n;

		Node_Link nl;
		nl.i = i;
		nl.in_link.clear();
		nl.out_link.clear();
		nl.nlink = 0;
        link_map[nid] = nl;
    }

	// setting link of map
	for (int64_t i=0; i<in_swc.row.size(); i++)
	{
        int64_t nid = in_swc.row.at(i).n;
		int64_t parent = in_swc.row.at(i).parent;

		if (parent>=0) // 0/1-based? valid link
		{
			Node_Link & nl = link_map[nid];
			nl.out_link.push_back(parent);
			nl.nlink ++;

			Node_Link & pnl = link_map[parent];
			pnl.in_link.push_back(nid);
            pnl.nlink ++;
        }else
        {
            bool single_node = true;
            for (int64_t j=0; j<in_swc.row.size(); j++)
            {
                if(in_swc.row.at(j).parent == nid)
                {
                    single_node = false;
                    break;
                }
            }
            if(single_node)
            {
                Node_Link & nl = link_map[nid];
                nl.out_link.push_back(nid);
                nl.nlink ++;
            }
        }
    }

	return link_map;
}

//091212 RZC
vector <V_NeuronSWC> decompose_V_NeuronSWC(V_NeuronSWC & in_swc)
{
	// map swc's index --> vector's index & in/out link
	Link_Map link_map = get_link_map(in_swc);

	vector <V_NeuronSWC> out_swc_segs;
	out_swc_segs.clear();

	// nchild as processed counter
	for (int64_t i=0; i<in_swc.row.size(); i++)
	{
		V_NeuronSWC_unit & cur_node = in_swc.row[i];
		Node_Link & nodelink = link_map[int64_t(cur_node.n)];
		cur_node.nchild = nodelink.nlink;
		//qDebug("#%d nlink = %d, in %d, out %d", int64_t(cur_node.n), nodelink.nlink, nodelink.in_link.size(), nodelink.out_link.size());
	}

	for (;;)
	{
		// check is all nodes processed
		int64_t n_removed = 0;
		for (int64_t i=0; i<in_swc.row.size(); i++)
		{
			V_NeuronSWC_unit & cur_node = in_swc.row[i];
			if (cur_node.nchild <=0) // is removed node ////// count down to 0
				n_removed ++;
		}
		if (n_removed >= in_swc.row.size())
		{
			//qDebug("split_V_NeuronSWC_segs is done.");
			break; //over, all nodes have been labeled to remove
		}

		// find a tip/out-branch/pure-out point as start point
		int64_t istart = -1;
		int64_t n_left = 0;
		int64_t i_left = -1;
		for (int64_t i=0; i<in_swc.row.size(); i++)
		{
			V_NeuronSWC_unit & cur_node = in_swc.row[i];
			Node_Link & nodelink = link_map[int64_t(cur_node.n)];

			if (cur_node.nchild <=0)
				continue; //skip removed point

			n_left++; //left valid point
			i_left = i;

			if ((nodelink.nlink ==1 && nodelink.in_link.size() == 0) // tip point (include single point)
			 || (nodelink.nlink >2 && nodelink.out_link.size() > 0) // out-branch point
			 || (nodelink.nlink ==2 && nodelink.in_link.size() == 0)) // pure-out point
			{
				istart = i;
				if (0)
                    cout << "start from " <<  int64_t(cur_node.n) << endl;
				break; //find a start point
			}
		}
		if (istart <0) //not find a start point
		{
			if (n_left)
				istart = i_left;
			else
				break;
		}

		// extract a simple segment
		V_NeuronSWC new_seg;
		new_seg.clear();
		//qDebug("decompose_V_NeuronSWC_segs: segment from node #%d", j);

		int64_t inext = istart;
		for (int64_t n=1; inext>=0; n++)
		{
			V_NeuronSWC_unit & cur_node = in_swc.row[inext];
			Node_Link & nodelink = link_map[int64_t(cur_node.n)];
			//qDebug("	link #%d", int64_t(cur_node.n));

			V_NeuronSWC_unit new_node = cur_node;
			new_node.n = n;
			new_node.parent = n+1; // link order as original order
			new_seg.row.push_back(new_node);

			if(cur_node.parent <0)    // root point ////////////////////////////
			{
				//qDebug("decompose_V_NeuronSWC_segs: segment end at root #%d", int64_t(cur_node.n));
				cur_node.nchild --;
				break; //over, a simple segment
			}
			else if (n>1 && nodelink.nlink >2)  // branch point (in link_map) ///////////
			{
				//qDebug("decompose_V_NeuronSWC_segs: segment end at branch #%d", int64_t(cur_node.n));
				cur_node.nchild --;
				break; //over, a simple segment
			}
			else if (n>1 && inext==istart)  // i_left point (a loop) ///////////
			{
				//qDebug("decompose_V_NeuronSWC_segs: segment end at branch #%d", int64_t(cur_node.n));
				cur_node.nchild --;
				break; //over, a simple segment
			}
			else  //(nodelink.nlink==2)   // path node ///////////////////////////////////////////////
			{
				//qDebug("decompose_V_NeuronSWC_segs: node #%d", j);
				int64_t parent = cur_node.parent;
				inext = link_map[parent].i; //// next point in seg
				cur_node.nchild = -1;  // label to remove
			}
		}

		if (new_seg.row.size()>0)//>=2)//? single point
		{
			new_seg.row[new_seg.row.size()-1].parent = -1; // set segment end
			char buf[10];
			new_seg.name = sprintf(buf,"%d", out_swc_segs.size()+1);
			new_seg.b_linegraph=true; //donot forget to do this //####################
			out_swc_segs.push_back(new_seg);
		}
	}

	return out_swc_segs;
}

//091027 RZC
bool reverse_V_NeuronSWC_inplace(V_NeuronSWC & in_swc)
{
	if (! in_swc.isLineGraph())
	{
		printf("Not a simple line-graph path! do not reverse.\n");
		return false;
	}

	V_NeuronSWC & out_swc = in_swc;

	bool order;
	for (int i=0; i<out_swc.row.size(); i++)
	{
		V_NeuronSWC_unit & node = out_swc.row[i];
		if (i==0)
		{
			order = (node.parent >0); // same link order as vector
		}

		if (order) // same --> revese
		{
			if (i==0)
			{
				node.parent = -1;
			}
			else
			{
				V_NeuronSWC_unit & prev_node = out_swc.row[i-1];
				node.parent = prev_node.n;
			}
		}
		else     // reverse --> same
		{
			if (i==out_swc.row.size()-1)
			{
				node.parent = -1;
			}
			else
			{
				V_NeuronSWC_unit & next_node = out_swc.row[i+1];
				node.parent = next_node.n;
			}
		}
	}
	return true;
}


bool verifyIsLineGraph(const V_NeuronSWC & in_swc) //to implement
{
//	vector <V_NeuronSWC_unit> * p_row = &(in_swc.row);
//	int64_t n=p_row->size(), i=0;
//	if (n<=0) return false;
//	int *h=0;
//	try{h=new int [n];}	catch(...)	return false;
//	for (i=0;i<n;i++) h[i]=0;
//
//	//compute the map table
//	for (i=0;i<n;i++)
//	{
//		h[i] = p_row->at.data[0];
//	}
//
//	//find the # parent array for all nodes
//
//		for (i=0;i<n;i++)
//		{
//			int64_t j=p_row->at(i).data[6];
//			if (j>)
//			if (>=0) h[p_row->at(i).data[6]]++;
//		}
//		if (h) {delete []h; h=0;}
//	}
	return false;
}

V_NeuronSWC_list split_V_NeuronSWC_simplepath(V_NeuronSWC & in_swc, int64_t nodeinseg_id)
{
	V_NeuronSWC_list out_swc_list;

	//first check if a line graph
	if (! in_swc.isLineGraph())
	{
		printf("Not a simple line-graph path! just copy instead of split.\n");
		out_swc_list.append(in_swc);
		return out_swc_list;
	}

	//now split
	//when the node is a tip, nothing should be done
	bool b_exist = (nodeinseg_id>0  &&  nodeinseg_id<in_swc.nrows()-1);

	if (!b_exist)  //if the node_id is invalid, then just copy
	{
		printf("invalid node index or tip of segment! just copy instead of split.\n");
		out_swc_list.append(in_swc);
		return out_swc_list;
	}

	if (b_exist)
	{
		V_NeuronSWC  res_swc1, res_swc2;
		int64_t j;

		for (j=0;j<=nodeinseg_id;j++)	res_swc1.append(in_swc.row.at(j)); //note to include the breaking point
		res_swc1.b_linegraph=true; //this must be set

		for (j=nodeinseg_id;j<in_swc.row.size();j++)	res_swc2.append(in_swc.row.at(j)); //note to include the breaking point
		res_swc2.b_linegraph=true; //this must be set

		// 091026 RZC: fixed the bug about link order
		if (in_swc.row[0].parent<1) // original root in first
		{
			res_swc2.row[0].parent = -1; //produce a root node for second
		}
		else // original root in second
		{
			res_swc1.row[nodeinseg_id].parent = -1; //produce a root node for first
		}

		out_swc_list.append(res_swc1);
		out_swc_list.append(res_swc2);
	}

	return out_swc_list;
}

V_NeuronSWC_list split_V_NeuronSWC_simplepath(V_NeuronSWC & in_swc, double xx, double yy, double zz)
{
	//now split
	int64_t nodeinseg_id = find_node_in_V_NeuronSWC(in_swc, xx,yy,zz);
	return split_V_NeuronSWC_simplepath(in_swc, nodeinseg_id);
}

int64_t find_node_in_V_NeuronSWC(V_NeuronSWC & in_swc, double x, double y, double z) //find the id of a node
{
	double eps = 1e-6;
	bool b_found=false;
	for (int64_t j=0;j<in_swc.row.size();j++)
	{
		const V_NeuronSWC_unit & node = in_swc.row.at(j);
		if (fabs(node.x-x)<eps
			&& fabs(node.y-y)<eps
			&& fabs(node.z-z)<eps
			)
		{
			b_found=true;
			return j;
		}
	}
	if (!b_found)
	{
		printf("can not find node contained the coordinate\n");
		return -1;
	}
}

int64_t find_seg_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, double x, double y, double z, int64_t & nodeinseg_id) //find the id of a seg
{
	double eps = 1e-6;
	bool b_found=false;
	for (int64_t k=0;k<swc_list.seg.size();k++)
	{
		int64_t nodeinseg_id = find_node_in_V_NeuronSWC(swc_list.seg[k], x,y,z);
		if (nodeinseg_id>=0)
		{
			b_found=true;
			return k;
		}
	}
	if (!b_found)
	{
		printf("can not find segment contained the coordinate\n");
		return -1;
	}
}

//int64_t find_seg_num_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t node_id) //find the id of a seg
//{
//	bool b_found=false;
////	for (int64_t k=0;k<swc_list.seg.size();k++)
////	{
////		for (int64_t j=0;j<swc_list.seg.at(k).row.size();j++)
////			if (swc_list.seg.at(k).row.at(j).data[0]==node_id)
////			{
////				b_found=true;
////				return k;
////			}
////	}
//
////	int64_t n=0;
////	for (int64_t k=0;k<swc_list.seg.size();k++)
////	{
////		if (node_id < n+swc_list.seg.at(k).row.size()) return k;
////		else n += swc_list.seg.at(k).row.size();
////	}
//	if (!b_found) return -1;
//}

bool change_type_in_seg_of_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id, int type)
{
	if (seg_id<0 || seg_id>= swc_list.seg.size())	return false;

	vector <V_NeuronSWC_unit> & row = (swc_list.seg[seg_id].row);
	for (int64_t j=0;j<row.size();j++)
	{
        row[j].type = type;
	}
	return true;
}

bool change_radius_in_seg_of_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id, double radius)
{
	if (seg_id<0 || seg_id>= swc_list.seg.size())	return false;

	vector <V_NeuronSWC_unit> & row = (swc_list.seg[seg_id].row);
	for (int64_t j=0;j<row.size();j++)
	{
		row[j].r = radius;
	}
	return true;
}

bool delete_seg_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id) //delete a seg in the V_NeuronSWC_list
{
	if (seg_id>=0 && seg_id<swc_list.seg.size()) 
	{
		swc_list.seg.erase(swc_list.seg.begin()+seg_id); 
		swc_list.last_seg_num=swc_list.seg.size(); //added on 2010-02-10. by PHC, guess probably is RZC forgot to add
		return true;
	}
	else return false;
}

double length_seg_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id) //find the length of a seg
{
	if (seg_id>=0 && seg_id<swc_list.seg.size())
	{
		return getLength_V_NeuronSWC(swc_list.seg.at(seg_id));
	}
	else return -1; //-1 is an invalid value as the dist must >=0
}

double getLength_V_NeuronSWC(V_NeuronSWC & subject_swc) //compute the length of a swc neuron
{
	map <int64_t, int64_t> sub_index_map = unique_V_NeuronSWC_nodeindex(subject_swc);
	double length=0;
	for (int64_t i=0;i<subject_swc.nrows();i++)
	{
		if (subject_swc.row.at(i).data[6]>=0)
		{
			int64_t parent_row = sub_index_map[subject_swc.row.at(i).data[6]];
			printf("now row = %d, parent row = %d \n", i, parent_row);
			double dx = subject_swc.row.at(i).data[2] - subject_swc.row.at(parent_row).data[2];
			double dy = subject_swc.row.at(i).data[3] - subject_swc.row.at(parent_row).data[3];
			double dz = subject_swc.row.at(i).data[4] - subject_swc.row.at(parent_row).data[4];
			length += sqrt(dx*dx+dy*dy+dz*dz);
		}
	}
	return length;
}


V_NeuronSWC join_V_NeuronSWC_vec(vector<V_NeuronSWC> & in_swc_vec)
{
	V_NeuronSWC out_swc;
	for (int64_t k=0;k<in_swc_vec.size();k++)
	{
		if (k==0) out_swc = in_swc_vec.at(0);
		else
		{
			if(!join_two_V_NeuronSWC(out_swc, in_swc_vec.at(k)))
			{
                printf("some error happened in join_V_NeuronSWC_vec() for seg [%d]\n",k);
			}
		}
	}
	out_swc.b_jointed = true; //091029 RZC
	return out_swc;
}

V_NeuronSWC join_segs_in_V_NeuronSWC_list(V_NeuronSWC_list & swc_list, int64_t seg_id_array[], int n_segs) //to implement
//merge several segs in V_NeuronSWC_list, indexed by seg_id_array (length is n_segs)
{
	V_NeuronSWC out_swc;
	if (!seg_id_array || n_segs<0 || n_segs>swc_list.nsegs())
	{
        printf("invalid parameters to join_segs_in_V_NeuronSWC_list()\n");
		return out_swc;
	}

	int64_t k;
	for (k=0;k<swc_list.nsegs();k++)
	{
		if (k==0) out_swc = swc_list.seg.at(k);
		else
        {
			if(!join_two_V_NeuronSWC(out_swc, swc_list.seg.at(k)))
			{
                printf("some error happened in join_segs_in_V_NeuronSWC_list() for seg [%d]\n",k);
			}
		}

		//save a temp swc file
//		char fname[128];
//		sprintf(fname, "/Users/pengh/temp/testaaa%d.swc", k);
//		FILE * fp = fopen(fname, "wt");
//		if (fp)
//		{
//			fprintf(fp, "#name %s\n", fname);
//			fprintf(fp, "#comment %s\n", "jointed neuron");
//
//			fprintf(fp, "#n,type,x,y,z,radius,parent\n");
//			for (int64_t i=0;i<out_swc.row.size();i++)
//			{
//				V_NeuronSWC_unit v = out_swc.row.at(i);
//				fprintf(fp, "%ld %ld %5.3f %5.3f %5.3f %5.3f %ld\n", int64_t(v.data[0]), int64_t(v.data[1]), v.data[2], v.data[3], v.data[4], v.data[5], int64_t(v.data[6]));
//			}
//
//			fclose(fp);
//		}


	}

	out_swc.b_jointed = true; //091029 RZC
	return out_swc;
}

bool join_two_V_NeuronSWC_old(V_NeuronSWC & destination_swc, V_NeuronSWC & subject_swc)
{
	if (subject_swc.nrows()<=0)
	{
        printf("subject_swc is empty in join_two_V_NeuronSWC()\n");
		return false;
	}

	if (destination_swc.nrows()<=0) {destination_swc = subject_swc; return true;}
	else
	{
		if (!simplify_V_NeuronSWC_nodeindex(destination_swc))
		{
            printf("inconsistent destination neuron structure found in join_two_V_NeuronSWC()\n");
			return false; //some error happens
		}

		V_NeuronSWC_unit v;
		int64_t i=0;
		destination_swc.printInfo();
		//one problem I don't consider now is that the nodes (3D locations) in a neuron reconstruction may duplicated. This should be simplified later too

		int64_t n0 = destination_swc.nrows();

		map <int64_t, int64_t> sub_index_map = unique_V_NeuronSWC_nodeindex(subject_swc);
		int64_t *ipos=0; ipos=new int64_t [subject_swc.nrows()];
		cout << "start ==================== map " << subject_swc.name.c_str() << endl;
		for (i=0;i<subject_swc.nrows();i++)
		{
			ipos[i] = destination_swc.getFirstIndexof3DPos(subject_swc.row.at(i));
			if (ipos[i]>=0) printf("**detect an overlapping node [sub ind = %d, dest_ind=%d]\n",i,ipos[i]);
			sub_index_map[subject_swc.row.at(i).data[0]] = (ipos[i]>=0) ? destination_swc.row.at(ipos[i]).data[0] : (sub_index_map[subject_swc.row.at(i).data[0]]+n0);
			cout << "map " <<subject_swc.row.at(i).data[0] << " to " << sub_index_map[subject_swc.row.at(i).data[0]] << endl;
		}
		cout <<"finish ==================== map" << endl;

		for (i=0;i<subject_swc.nrows();i++)
		{
			v = subject_swc.row.at(i);
			cout << i << "before : " << v.data[0] << " " << v.data[1] << " "
                 << v.data[2] << " " << v.data[3] << " " << v.data[4] <<" "
                 << v.data[5] << " " << v.data[6] << endl;

			v.data[0] = sub_index_map[v.data[0]];
			if (ipos[i]>=0) //overlap
			{
				if (v.data[6]>=0) //non-root
				{
					if (ipos[int64_t(v.data[6])]<0) //if the parent overlap as well, then do no add
					{
						v.data[6] = sub_index_map[v.data[6]]; //if the overalp node is not a root (thus has its own parent, then inherit)
						destination_swc.append(v); //in this case, still need to add this swc unit into the dest_swc
					}
				}
				else //root
				{
					v.data[6] = destination_swc.row.at(ipos[i]).data[6]; //if the overlap node is a root for subject seg, then copy the parent in dest_swc
					cout << "**omit one overlap node **" << endl;
					//in this case, no need to add this swc unit into the dest_swc
				}
			}
			else //if not an overlap node, should just add (but the node indexes should be updated)
			{
				if (v.data[6]>=0)
				{
					v.data[6] = sub_index_map[v.data[6]]; //if the overalp node is not a root (thus has its own parent, then inherit)
					destination_swc.append(v);
				}
				else
				{
					//if the overlap node is a root for subject seg, then copy the parent in dest_swc
					cout << "**copy a non-overlapping root node **" << endl;
					destination_swc.append(v);	//in this case, still need to add as it is not overlap
				}
			}
			cout << i << "after : " << v.data[0] << " " << v.data[1] << " "
                 << v.data[2] << " " << v.data[3] << " " << v.data[4] << " "
                 << v.data[5] << " " << v.data[6] << endl;
		}

		//search for all root nodes
		destination_swc.printInfo();
		vector < vector<int64_t> > des_rootnode_overlap_info;
		vector <int64_t> overlap_root_array;
		V_NeuronSWC newswc = destination_swc; newswc.row.clear();
		for (i=0; i<destination_swc.row.size();i++) //note that I search the entire merged destination, so that the part originally from subject will get searched again
		{
			v = destination_swc.row.at(i);
			if (v.parent<0)
			{
				vector <int64_t> mypos = destination_swc.getAllIndexof3DPos(v, i);
				if (mypos.size()>1)
				{
					printf("Detect a root node that has the same coordinates of more than one other neuron nodes. Your data may have error. \n");
				}

				if (mypos.size()>0) //overlap root, remove
				{
					vector <int64_t> mypos1;
					for (int64_t j=0;j<mypos.size();j++)
					{
						if (destination_swc.row.at(mypos.at(j)).parent<0) //if it is a root
						{
							if (destination_swc.row.at(mypos.at(j)).n!=v.n) //then it is not a redundant row, -- maybe there is some other nodes link to this one. then do not add it to the newswc
							{
								destination_swc.row.at(mypos.at(j)).n = v.n;
							}
						}
						else //the overlap one is not a root
						{
							if (destination_swc.row.at(mypos.at(j)).n!=v.n)
								mypos1.push_back(destination_swc.row.at(mypos.at(j)).n);
							//if equal, do nothing (which means remove this root-row. In this case, no node mapping is needed, but the root should be deleted
						}
					}

					if (mypos1.size()>0)
					{
						overlap_root_array.push_back(v.n);
						des_rootnode_overlap_info.push_back(mypos1);
						//remove the current row by not copy to the new structure
						cout << "remove overlap root node" << v.n << endl;
					}
					//newswc.append(v);
				}
				else //non-overlapping root, keep
				{
					newswc.append(v);
				}
			}
			else //non_root, keep
			{
				newswc.append(v);
			}
		}

		for (i=0;i<overlap_root_array.size();i++)
            cout << "overlapping root array " << overlap_root_array.at(i) << endl;

		for (i=0;i<newswc.row.size();i++) //now merge & update
		{
			v = newswc.row.at(i);
			int64_t cur_parent = int64_t(v.parent);
			int64_t mytmp = value_in_vector(overlap_root_array, cur_parent);
			if (mytmp>=0)
			{
                cout << "cur node " << v.n << " cur parent " << v.parent
                     << "new parent " << des_rootnode_overlap_info.at(mytmp).at(0) << endl;
				newswc.row.at(i).parent = des_rootnode_overlap_info.at(mytmp).at(0);

				for (int64_t j=1;j<des_rootnode_overlap_info.at(mytmp).size();j++)
				{
					v.parent = des_rootnode_overlap_info.at(mytmp).at(j);
					newswc.append(v); //append some new records if one root mapped to multiple other non-root nodes
				}
			}
		}
		destination_swc = newswc; //update

		//save a temp swc file
		FILE * fp = fopen("/Users/pengh/temp/testaaa.swc", "wt");
		if (fp)
		{
			fprintf(fp, "#name %s\n", "test123.swc");
			fprintf(fp, "#comment %s\n", "jointed neuron");

			fprintf(fp, "#n,type,x,y,z,radius,parent\n");
			for (i=0;i<destination_swc.row.size();i++)
			{
				v = destination_swc.row.at(i);
				fprintf(fp, "%ld %ld %5.3f %5.3f %5.3f %5.3f %ld\n", int64_t(v.data[0]), int64_t(v.data[1]), v.data[2], v.data[3], v.data[4], v.data[5], int64_t(v.data[6]));
			}

            fclose(fp);
		}

		if(ipos) {delete []ipos; ipos=0;}
	}

	destination_swc.b_jointed = true; //091029 RZC
	return true;
}


bool join_two_V_NeuronSWC(V_NeuronSWC & destination_swc, V_NeuronSWC & subject_swc)
{
	if (subject_swc.nrows()<=0) {printf("subject_swc is empty in join_two_V_NeuronSWC()\n"); return false;}
	if (destination_swc.nrows()<=0) {destination_swc = subject_swc; return true;}

	//step 0: to simplify the des neuron, to make the joining simpler
	if (!simplify_V_NeuronSWC_nodeindex(destination_swc))
	{
		printf("inconsistent destination neuron structure found in join_two_V_NeuronSWC()\n");
		return false; //some error happens
	}

    //first just concatenate records
	V_NeuronSWC_unit v, vp;
	int64_t i=0, j=0;

	int64_t nr_des = destination_swc.nrows();
	for (i=0;i<subject_swc.nrows();i++)
	{
        v = subject_swc.row.at(i);
		v.n += nr_des;
		if (v.parent>=0) v.parent += nr_des;
		destination_swc.append(v);
	}
    nr_des = destination_swc.nrows();

	//then produce the adjacency matrix
	vector<V_NeuronSWC_coord> unpos = destination_swc.unique_ncoord();
	int64_t N = unpos.size();

	map <int64_t, int64_t> nid_row_lut;
	for (i=0; i<nr_des; i++) { nid_row_lut[destination_swc.row.at(i).n] = i; }

	map <int64_t, int64_t> row_unid_lut;
	for (i=0; i<nr_des; i++)
	{
		V_NeuronSWC_coord c = destination_swc.row.at(i).get_coord();
		int64_t ic = value_in_vector(unpos, c);
		if (ic<0)
		{
			printf("Indexing error! Check data in join_two_V_NeuronSWC.\n");
			return false;
		}
		row_unid_lut[i] = ic;
	}

	vector < vector <int64_t> > adjm;
	vector <int64_t> plist;
	for (i=0; i<N; i++) adjm.push_back(plist); //initialize as empty

	vector < vector <float> > adjm_radius;
	vector <float> tmpr;
	for (i=0; i<N; i++) adjm_radius.push_back(tmpr); //initialize as empty

    vector < vector <float> > adjm_type;
    vector <float> tmpt;
	for (i=0; i<N; i++) adjm_type.push_back(tmpt); //initialize as empty

	for (i=0; i<nr_des; i++)
	{
        v = destination_swc.row.at(i);
        int64_t iv = row_unid_lut[i];

		if (v.parent>=0)
		{
			vp = destination_swc.row.at(nid_row_lut[v.parent]);
			int64_t ip =  row_unid_lut[nid_row_lut[v.parent]];
			adjm.at(iv).push_back(ip); //add it anyway, will find unique one later. Indeed even unnecessary to find the unique ones, because the repeated record will merge automatically
			adjm_radius.at(iv).push_back(float(v.r));
            adjm_type.at(iv).push_back(float(v.type));
			//need to record neuron edge type and radius as well later
		}
		else
		{
			//do nothing about the parent, i.e. this is a root for this record, and thus do not add to the parent list at all
			//but still update type and radius
			adjm_radius.at(iv).push_back(float(v.r));
            adjm_type.at(iv).push_back(float(v.type));
        }
	}

	//no need to find unique ones

	//then re-produce a new swc
	V_NeuronSWC newswc = destination_swc; newswc.row.clear();
	for (i=0;i<N;i++)
	{
		vector <int64_t> & cur_pa = adjm.at(i);
		vector <float> & cur_r = adjm_radius.at(i);
        vector <float> & cur_t = adjm_type.at(i);
		V_NeuronSWC_coord & cur_coord = unpos.at(i);

		if (cur_pa.size()<=0) //a root
		{
			v.n = i;
            v.type = double(cur_t.back()); //direct use the last one, which means the later one overwrite the first one
			v.x = cur_coord.x;
			v.y = cur_coord.y;
			v.z = cur_coord.z;
			v.r = cur_r.back(); //direct use the last one, which means the later one overwrite the first one
			v.parent = -1;
			newswc.append(v);
		}
		else
		{
			for (j=0;j<cur_pa.size();j++)
			{
				v.n = i;
                v.type = double(cur_t.at(j)); //allow he later one overwrite the earlier one(s) if there are redundant records
				v.x = cur_coord.x;
				v.y = cur_coord.y;
				v.z = cur_coord.z;
				v.r = cur_r.at(j); //so I allow the later one overwrite the earlier one(s) if there are redundant records
				v.parent = cur_pa.at(j);
				newswc.append(v);
			}
        }
	}


	//copy the results
	destination_swc = newswc;
    destination_swc.b_jointed = true;
	return true;
}




map <int64_t,int64_t> unique_V_NeuronSWC_nodeindex(V_NeuronSWC & in_swc)
{
	map <int64_t, int64_t> index_map;
	const int64_t n = in_swc.nrows();
	for (int64_t i=0;i<n;i++)
	{
		index_map[int64_t(in_swc.row.at(i).data[0])]=i;
	}
	return index_map;
}

bool simplify_V_NeuronSWC_nodeindex(V_NeuronSWC & my_swc)
 //map the node index of a swc neuron to the range 1~N (N is the number of unique indexes)
{
	map <int64_t, int64_t> index_map = unique_V_NeuronSWC_nodeindex(my_swc);
	map <int64_t, int64_t>::iterator it;
	V_NeuronSWC out_swc = my_swc;
	for (int64_t i=0; i< my_swc.nrows(); i++)
	{
		out_swc.row.at(i).n = index_map[my_swc.row.at(i).n];
		if (my_swc.row.at(i).parent>=0)
		{
			it = index_map.find(my_swc.row.at(i).parent);
			if (it==index_map.end())
			{
				cout << i << " " << my_swc.row.at(i).parent << endl;
                cout << "detected an inconsistent case in simplify_V_NeuronSWC_nodeindex()" << endl;
				return false; //this indicates the SWC structure has inconsistency
			}
			else out_swc.row.at(i).parent = (*it).second;
		}
	}
	my_swc = out_swc; //now update. I use this slower version to ensure the my_swc will NOT be partially damaged even when its content is not consistent
	return true;
}


V_SWCNodes V_NeuronSWC::unique_nodes_info()
{
	V_SWCNodes res;	if (row.size()<1) return res;
	res.nid_array.push_back(row.at(0).n);
	res.ncoord_array.push_back(row.at(0).get_coord());
	res.nid_ipos_lut[res.nid_array.at(0)] = 0;

	for (int64_t i=1;i<row.size();i++)
	{
		int64_t cur_nid = int64_t(row.at(i).n);
		V_NeuronSWC_coord cur_ncoord = row.at(i).get_coord();

		int64_t ind_nid = value_in_vector(res.nid_array, cur_nid);
		int64_t ind_nccord = value_in_vector(res.ncoord_array, cur_ncoord);
		if (ind_nid==ind_nccord && ind_nid>=0 && ind_nccord>=0)
		{
			res.nid_array.push_back(cur_nid);
			res.ncoord_array.push_back(cur_ncoord);
			res.nid_ipos_lut[cur_nid] = ind_nid;
		}
	}
	return res;
}



vector<int64_t> V_NeuronSWC::unique_nid()
{
	vector<int64_t> res;	if (row.size()<1) return res;
	res.push_back(row.at(0).n);
	for (int64_t i=1;i<row.size();i++)
	{
		int64_t cur_nid = int64_t(row.at(i).n);
		if (value_in_vector(res, cur_nid)<0)
			res.push_back(cur_nid);
	}
	return res;
}

int64_t V_NeuronSWC::n_unique_nid()
{
	vector<int64_t> res = unique_nid();
	return res.size();
}

vector<V_NeuronSWC_coord> V_NeuronSWC::unique_ncoord()
{
	vector<V_NeuronSWC_coord> res;	if (row.size()<1) return res;
	V_NeuronSWC_coord cur_coord; cur_coord.set(row.at(0).x, row.at(0).y, row.at(0).z);
	res.push_back(cur_coord);
	for (int64_t i=1;i<row.size();i++)
	{
		cur_coord.set(row.at(i).x, row.at(i).y, row.at(i).z);
		if (value_in_vector(res, cur_coord)<0)
			res.push_back(cur_coord);
	}
	return res;
}

int64_t V_NeuronSWC::n_unique_ncoord()
{
	vector<V_NeuronSWC_coord> res = unique_ncoord();
	return res.size();
}




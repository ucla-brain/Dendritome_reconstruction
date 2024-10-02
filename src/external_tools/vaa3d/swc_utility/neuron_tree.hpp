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

//basic_surf_objs.h
//defining the basic types of surface objects
//by Hanchuan Peng. 2009-03-06. This is the first step to unify all basic surface object data types in different modules
//
//090605: merge with p_objectfile.h
//090706: add neuronswc and marker write functions here
//130103: add a few constructor functions for ImageMarker

// modified from vaa3d basic_surf_objs.h  mzhu 05/23/2019

#ifndef MCP3D_VAA3D_NEURON_TREE_HPP
#define MCP3D_VAA3D_NEURON_TREE_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include "xyz.hpp"

// .swc neurons and other graph-describing files
// 5/30/2019 change xyzr type from float to double mzhu. in MyMarker these values
// are double type. swc files written from MyMarker and read to NeuronSWC
// may have floating point out of range issue
struct NeuronSWC
{
    int64_t n;
	int type;			// 0-Undefined, 1-Soma, 2-Axon, 3-Dendrite, 4-Apical_dendrite, 5-Fork_point, 6-End_point, 7-Custom
	double x, y, z;		// point coordinates
    
    union{
    double r;			// radius
    double radius;
    };
    
    union{
	int64_t pn;				// previous point index (-1 for the first point)
	int64_t parent;				// previous point index (-1 for the first point)
    };
    
    int64_t level; //20120217, by PHC. for ESWC format
    std::vector<float> fea_val; //20120217, by PHC. for ESWC format

	int64_t seg_id; //this is reused for ESWC format, 20120217, by PHC
    int64_t nodeinseg_id; //090925, 091027: for segment editing

	NeuronSWC () {type = 0; n = pn = 0; x=y=z=r=0; seg_id=-1; nodeinseg_id=0; fea_val=std::vector<float>{}; level=-1;}
};

// .neuron trees

struct NeuronTree
{
	std::vector<NeuronSWC> listNeuron;
	std::unordered_map<int64_t, int64_t>  hashNeuron;
	std::string file_name;

    int64_t n;

    NeuronTree()
    {
        n = 0;
        listNeuron.clear(); 
        hashNeuron.clear(); 
        file_name="";
    }

    void deepCopy(const NeuronTree p)
    {
        n=p.n;

        file_name = p.file_name;
        hashNeuron.clear();

        for (int i =0; i< p.listNeuron.size();i++){
            NeuronSWC S;
            S.n = p.listNeuron[i].n;
            S.type = p.listNeuron[i].type;
            S.x = p.listNeuron[i].x;
            S.y= p.listNeuron[i].y;
            S.z = p.listNeuron[i].z;
            S.r = p.listNeuron[i].r;
            S.pn = p.listNeuron[i].pn;
            S.seg_id = p.listNeuron[i].seg_id;
            S.fea_val = p.listNeuron[i].fea_val;
            listNeuron.push_back(S);
            hashNeuron.insert(std::make_pair(S.n, listNeuron.size() - 1));
        }
    }

	void copy(const NeuronTree & p)
	{
		n=p.n;
		listNeuron = p.listNeuron;
		hashNeuron = p.hashNeuron;
		file_name     = p.file_name;
	}
	void copyGeometry(const NeuronTree & p)
	{
		if (p.listNeuron.size()!=listNeuron.size()) return;

		NeuronSWC *p_tmp;
		for (size_t i=0; i < listNeuron.size(); i++)
		{
			p_tmp = (NeuronSWC *)(&(listNeuron.at(i)));
			//qDebug()<<"before:"<<p_tmp->x<<p_tmp->y<<p_tmp->z<<p_tmp->r;
			p_tmp->x = p.listNeuron.at(i).x;
			p_tmp->y = p.listNeuron.at(i).y;
			p_tmp->z = p.listNeuron.at(i).z;
			p_tmp->r = p.listNeuron.at(i).r;
			//qDebug()<<"src:"<<p.listNeuron.at(i).x<<p.listNeuron.at(i).y<<p.listNeuron.at(i).z<<p.listNeuron.at(i).r;
			//qDebug()<<"after:"<<p_tmp->x<<p_tmp->y<<p_tmp->z<<p_tmp->r;
		}
	}
	bool projection(int axiscode=3) //axiscode, 1 -- x, 2 -- y, 3 -- z, 4 -- r
	{
		if (axiscode!=1 && axiscode!=2 && axiscode!=3 && axiscode!=4)
            return false;
		NeuronSWC *p_tmp;
		for (int64_t i=0;i<listNeuron.size();i++)
		{
			p_tmp = (NeuronSWC *)(&(listNeuron.at(i)));
			//qDebug()<<"before:"<<p_tmp->x<<p_tmp->y<<p_tmp->z<<p_tmp->r;
			if (axiscode==1) p_tmp->x = 0;
			else if (axiscode==2) p_tmp->y = 0;
			else if (axiscode==3) p_tmp->z = 0;
			else if (axiscode==4) p_tmp->r = 0.5;
		}
		return true;
	}
};

//general operators
inline bool operator==(NeuronTree& a, NeuronTree& b)
{
	return a.file_name == b.file_name;
}


#endif // MCP3D_VAA3D_SD_NEURON_TREE_HPP


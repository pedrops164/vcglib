/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/*! \file trimesh_base.cpp
\ingroup code_sample

\brief the minimal example of using the lib

This file contain a minimal example of the library, showing how to load a mesh and how to compute per vertex normals on it.

*/

#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>
#include "my_mesh.h"
#include "main.h"
#include <iostream>
#include <iomanip>

using namespace vcg;
using namespace std;

int main( int argc, char **argv )
{
    if(argc<2)
    {
      printf("Usage ./prog meshfilename.off\n");
      return -1;
    }
    MyMesh m;
    
    cout << std::fixed;
    cout << std::setprecision(2);
    
    if(tri::io::ImporterOFF<MyMesh>::Open(m,argv[1])!=tri::io::ImporterOFF<MyMesh>::NoError)
    {
      printf("Error reading file  %s\n",argv[1]);
      exit(0);
    }
    
    tri::RequirePerVertexNormal(m);
    tri::UpdateTopology<MyMesh>::FaceFace(m);
    HoleFill<MyMesh>::my_hole_fill(m, "hole_filled.off");
    return 0;
}


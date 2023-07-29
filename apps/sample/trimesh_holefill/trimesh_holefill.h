#include "my_mesh.h"

using namespace vcg;
using namespace std;

#include "triangulation_alg.h"
#include "refinement_alg.h"
#include "fairing_alg.h"
#include "alg_util.h"
#include "map_EF.h"

template<class MeshType>
class HoleFill {

public:
	typedef typename MeshType::VertexType VertexType;
	typedef typename MeshType::FaceType FaceType;
	static void my_hole_fill(MeshType& input_mesh, const char* output_mesh) {
		tri::UpdateFlags<MeshType>::VertexClearV(input_mesh);
		int border_loop_count = 0;
		int border_vertex_count = 0;
		vector<MeshType*> patched_holes;
		// search for a boundary face
		for (auto fi = input_mesh.face.begin(); fi != input_mesh.face.end(); ++fi) {
			for (int i = 0; i < 3; ++i) {
				if (face::IsBorder(*fi, i) && !(*fi).V(i)->IsV()) {// mark visited the vertices of the border edge
					vector<VertexType*> borderVertices;
					vector<FaceType*> borderFaces; //vector of faces on the border
					EdgeFaceMap<MeshType> map;
					(*fi).V(i)->SetV();
					// now start a walk along the border
					borderVertices.push_back((*fi).V0(i));
					borderFaces.push_back(&*fi);
					border_vertex_count += 1;
					//++border_loop_count;
					face::Pos<FaceType> start(&*fi, i);
					face::Pos<FaceType> cur = start;
					do
					{
						assert(cur.IsBorder());
						do
						{
							cur.FlipE();
							cur.FlipF();
						} while (!cur.IsBorder());

						VertexType* previous = cur.V();
						cur.FlipV();
						VertexType* current = cur.V();
						map.add(previous, current, &*cur.F());
						if (!cur.V()->IsV()) {
							cur.V()->SetV();
							border_vertex_count++;
							borderVertices.push_back(cur.V());
							borderFaces.push_back(&*cur.F());
						}
						cur.F()->C() = Color4b::Red;
					} while (cur != start);
					//here one hole has been traversed all the way
					MeshType* final_patch = nullptr;
					MeshType* cpy = nullptr;

					//Triangulation Algorithm
					cout << "before triangulation" << endl;
					TriangulationAlg<MeshType> ta(borderVertices, map);
					final_patch = &ta.algorithm();
					cpy = AlgUtil<MeshType>::mesh_copy(input_mesh);
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, *final_patch);
					tri::io::ExporterOFF<MeshType>::Save(*cpy, ("hole" + to_string(border_loop_count) + "triangulation.off").c_str(), tri::io::Mask::IOM_FACECOLOR);
					AlgUtil<MeshType>::free_mesh(*cpy);


					//Refining Algorithm
					cout << "before refining" << endl;
					RefinementAlg<MeshType> ra(*final_patch, borderVertices, borderFaces);
					final_patch = &ra.algorithm();
					cpy = AlgUtil<MeshType>::mesh_copy(input_mesh);
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, *final_patch);
					tri::io::ExporterOFF<MeshType>::Save(*cpy, ("hole" + to_string(border_loop_count) + "refined.off").c_str(), tri::io::Mask::IOM_FACECOLOR);
					AlgUtil<MeshType>::free_mesh(*cpy);

					//Fairing Algorithm
					cout << "before fairing" << endl;
					FairingAlg<MeshType> fa(*final_patch, borderVertices, borderFaces);
					final_patch = &fa.algorithm();
					cpy = AlgUtil<MeshType>::mesh_copy(input_mesh);
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, *final_patch);
					tri::io::ExporterOFF<MeshType>::Save(*cpy, ("hole" + to_string(border_loop_count) + "faired.off").c_str(), tri::io::Mask::IOM_FACECOLOR);
					AlgUtil<MeshType>::free_mesh(*cpy);


					patched_holes.push_back(final_patch);
					border_loop_count++;
				}
			}
		}
		printf("Found %i border loops\n", border_loop_count);

			   //new vertices indexes correspond with the indexes of the border vertices
		for (MeshType* patched_hole : patched_holes) {
			vcg::tri::Append<MeshType, MeshType>::Mesh(input_mesh, *patched_hole);
		}
		// save the mesh with the border faces colored in red
		//tri::io::ExporterOFF<MeshType>::Save(input_mesh, output_mesh, tri::io::Mask::IOM_FACECOLOR);
	}

};

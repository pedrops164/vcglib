
#include "triangulation_alg.h"
#include "refinement_alg.h"
#include "fairing_alg.h"

template<class MeshType>
class HoleFill {

public:
	typedef typename MeshType::VertexType VertexType;
	typedef typename MeshType::FaceType FaceType;
	typedef typename vcg::face::Pos<FaceType> PosType;
	static void my_hole_fill(MeshType& surrounding_mesh, MeshType& patching_mesh, const char* output_mesh) {
		vcg::tri::UpdateFlags<MeshType>::VertexClearV(surrounding_mesh);
		int border_loop_count = 0;
		std::vector<MeshType*> patched_holes;
		// search for a boundary face
		auto begin = surrounding_mesh.face.begin();
		auto end = surrounding_mesh.face.end();
		for (auto fi = begin; fi != end; ++fi) {
			for (int i = 0; i < 3; ++i) {
				if (vcg::face::IsBorder(*fi, i) && !(*fi).V(i)->IsV()) {// mark visited the vertices of the border edge
					std::vector<VertexType*> borderVertices;
					std::vector<FaceType*> borderFaces; //vector of faces on the border
					(*fi).V(i)->SetV();
					// now start a walk along the border
					borderVertices.push_back((*fi).V0(i));
					borderFaces.push_back(&*fi);
					//++border_loop_count;
					vcg::face::Pos<FaceType> start(&*fi, i);
					vcg::face::Pos<FaceType> cur = start;
					std::vector<PosType> vv;
					do
					{
						assert(cur.IsBorder());
						vv.push_back(cur);
						do
						{
							cur.FlipE();
							cur.FlipF();
						} while (!cur.IsBorder());

						VertexType* previous = cur.V();
						cur.FlipV();
						VertexType* current = cur.V();
						if (!cur.V()->IsV()) {
							cur.V()->SetV();
							borderVertices.push_back(cur.V());
							borderFaces.push_back(&*cur.F());
						}
						cur.F()->C() = vcg::Color4b::Red;
					} while (cur != start);
					//here one hole has been traversed all the way
					patching_mesh.Clear();
					MeshType* cpy = nullptr;

					//Triangulation Algorithm
					std::cout << "before triangulation" << std::endl;
					TriangulationAlg<MeshType> ta(patching_mesh, vv);
					ta.algorithm();
					cpy = new MeshType();
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, surrounding_mesh);
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, patching_mesh);
					vcg::tri::io::ExporterOFF<MeshType>::Save(*cpy, ("hole" + std::to_string(border_loop_count) + "triangulation.off").c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);


					//Refining Algorithm
					std::cout << "before refining" << std::endl;
					RefinementAlg<MeshType> ra(patching_mesh, borderVertices, borderFaces);
					ra.algorithm();
					cpy = new MeshType();
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, surrounding_mesh);
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, patching_mesh);
					vcg::tri::io::ExporterOFF<MeshType>::Save(*cpy, ("hole" + std::to_string(border_loop_count) + "refined.off").c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);

					//Fairing Algorithm
					std::cout << "before fairing" << std::endl;
					FairingAlg<MeshType> fa(patching_mesh, borderVertices, borderFaces);
					fa.algorithm();
					cpy = new MeshType();
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, surrounding_mesh);
					vcg::tri::Append<MeshType, MeshType>::Mesh(*cpy, patching_mesh);
					vcg::tri::io::ExporterOFF<MeshType>::Save(*cpy, ("hole" + std::to_string(border_loop_count) + "faired.off").c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);


					//patched_holes.push_back(patched_mesh);
					vcg::tri::Append<MeshType, MeshType>::Mesh(surrounding_mesh, patching_mesh);
					border_loop_count++;
				}
			}
		}
		printf("Found %i border loops\n", border_loop_count);

			   //new vertices indexes correspond with the indexes of the border vertices
		//for (MeshType* patched_hole : patched_holes) {
		//	vcg::tri::Append<MeshType, MeshType>::Mesh(surrounding_mesh, *patched_hole);
		//}
		// save the mesh with the border faces colored in red
		//vcg::tri::io::ExporterOFF<MeshType>::Save(surrounding_mesh, output_mesh, vcg::tri::io::Mask::IOM_FACECOLOR);
	}

};

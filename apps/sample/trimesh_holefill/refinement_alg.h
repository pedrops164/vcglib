#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <vcg/complex/algorithms/update/color.h>

template<class MeshType>
class RefinementAlg {
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
private:
	double alpha;
	MeshType& mesh;
	std::vector<VertexType*> bV; //border vertices
	std::vector<FaceType*> bF; //border faces
	std::vector<double> scale; //vector of scale for each vertex. scale[i] is the scale of the i_th vertexs

public:

	RefinementAlg(MeshType& m, std::vector<VertexType*> bV, std::vector<FaceType*> bF): mesh(m), bV(bV), bF(bF) {
		alpha = sqrt(2);

		vcg::tri::RequirePerVertexNormal(mesh); //needed?

		size_t n = bV.size();		
		vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
	}

	void algorithm() {

			   //step 1
		/*  For each vertex v_i on the hole boundary, compute the
			scale attribute sigma(v_i) as the average length of the
			edges that are adjacent to v_i in the surrounding mesh.
			Initialize the patching mesh as the given hole
			triangulation.
		*/
		for (int i = 0; i < bV.size(); i++) {
			//iterate over all the faces on the border
			vcg::face::Pos<FaceType> cur(bF[i], bV[i]);
			if (!cur.IsBorder()) {
				//if it's not a border, switch the edge (to go to the border)
				cur.FlipE();
				assert(cur.IsBorder());
			}
			vcg::Point3f orig = cur.V()->P();
			double sum_distance = 0;
			int count = 0;
			do {
				vcg::Point3f dist = cur.VFlip()->P();
				sum_distance += vcg::Distance(orig, dist);
				count++;
				cur.FlipE();
				cur.FlipF();

			} while (!cur.IsBorder());
			//calculate last edge length
			vcg::Point3f dest = cur.VFlip()->P();
			sum_distance += vcg::Distance(orig, dest);
			count++;

			double avg_distance = sum_distance / count;
			scale.push_back(avg_distance);
			//every vertex on the border is mapped to its scale attribute
		}

		//step 2
		/*
		For each triangle in the patching mesh,
		compute the centroid v_c and the corresponding scale
		attribute. For m=i,j,k, if the alpha condition is met, then replace the triangle
		with the three smaller triangles in the patching mesh, and relax the edges of the original triangle.
		If the outer edges are on the border of the patching mesh, these are not relaxed because
		it would mean that the faces on the surrounding mesh would also be modified.

		 returns true if triangles were added , false otherwise
		 */
		auto step2 = [&]() {
			bool triangles_added = false;
		
			std::vector<int> faces_indices;
			for (auto fi = mesh.face.begin(); fi != mesh.face.end(); fi++) if (!fi->IsD()) {
				VertexType* v0 = fi->V(0);
				VertexType* v1 = fi->V(1);
				VertexType* v2 = fi->V(2);
				vcg::Point3f v0_p = v0->P();
				vcg::Point3f v1_p = v1->P();
				vcg::Point3f v2_p = v2->P();
				int v0_index = vcg::tri::Index(mesh, v0);
				int v1_index = vcg::tri::Index(mesh, v1);
				int v2_index = vcg::tri::Index(mesh, v2);
		
				vcg::Point3f centroid_cords = (v0_p + v1_p + v2_p) / 3;
				double scale_v0 = scale[v0_index];
				double scale_v1 = scale[v1_index];
				double scale_v2 = scale[v2_index];
				double scale_centroid = (scale_v0 + scale_v1 + scale_v2) / 3;
		
				//definition of alpha condition
				auto check = [&](VertexType* v) {
					vcg::Point3f dest = v->P();
					double dist = vcg::Distance(centroid_cords, dest);
					double value = alpha * dist;
					if (value > scale_centroid && value > scale[vcg::tri::Index(mesh, v)]) return true;
					return false;
				};
		
		
				//For m=i,j,k, if the alpha condition is met, subdivide triangle
				if (check(v0) || check(v1) || check(v2)) {
					triangles_added = true;
					FaceType* f = &*fi;
					faces_indices.push_back(vcg::tri::Index(mesh, f));
					scale.push_back(scale_centroid);
				}
				else {
					//otherwise we don't subdivide the current triangle
				}
			}
		
			//for each face to be split, we add three faces to the mesh, and we delete the current one
			auto f_p = vcg::tri::Allocator<MeshType>::AddFaces(mesh, faces_indices.size() * 3);
		
			//for each face to be added, one vertex is added (the centroid of each face)
			auto v_p = vcg::tri::Allocator<MeshType>::AddVertices(mesh, faces_indices.size());
		
			for (int face_idx : faces_indices) {
				FaceType* fi = &(mesh.face[face_idx]);
				vcg::tri::Allocator<MeshType>::DeleteFace(mesh, *fi); //deletes old face
				
				VertexType* v0 = fi->V(0);
				VertexType* v1 = fi->V(1);
				VertexType* v2 = fi->V(2);
				vcg::Point3f v0_p = v0->P();
				vcg::Point3f v1_p = v1->P();
				vcg::Point3f v2_p = v2->P();
				int v0_index = vcg::tri::Index(mesh, v0);
				int v1_index = vcg::tri::Index(mesh, v1);
				int v2_index = vcg::tri::Index(mesh, v2);
		
				vcg::Point3f centroid_cords = (v0_p + v1_p + v2_p) / 3;
				VertexType* centroid = &*v_p;
				v_p++;
				centroid->P() = centroid_cords;
		
				FaceType* new_face0 = &*f_p;
				f_p++;
				FaceType* new_face1 = &*f_p;
				f_p++;
				FaceType* new_face2 = &*f_p;
				f_p++;
				new_face0->V(0) = centroid;
				new_face0->V(1) = v2;
				new_face0->V(2) = v0;
				new_face0->C() = vcg::Color4b::White;

				new_face1->V(0) = centroid;
				new_face1->V(1) = v1;
				new_face1->V(2) = v2;
				new_face1->C() = vcg::Color4b::White;

				new_face2->V(0) = centroid;
				new_face2->V(1) = v0;
				new_face2->V(2) = v1;
				new_face2->C() = vcg::Color4b::White;
			}
			vcg::tri::Allocator<MeshType>::CompactFaceVector(mesh);
		
			//after all vertices have been added to new mesh, we update the topology
			vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
			//new_mesh = aux_mesh;
			return triangles_added;
		};

		//step 4
		/*
		Relax all interior edges of the patching mesh.
		*/
		auto step4 = [&]() {
			bool relaxed_edge = false;

			/*
			 * This lambda function receives a half edge position, and tries to relax the edge. It returns true if the edge
			 * was relaxed, false otherwise. If the half edge is a border, it just returns false.
			 */
			auto relax_edge = [&](vcg::face::Pos<FaceType> hedge) {
				if (hedge.IsBorder()) {
					return false; //we don't relax the edge, so we return false
				}

				//Here we calculate the 4 vertices involved and their coordinates
				int edge_to_swap = hedge.E();
				auto aux1 = hedge;
				aux1.FlipE();
				aux1.FlipV();
				FaceType* current_face = aux1.F();
				VertexType* v1 = aux1.V(); //internal vertex
				vcg::Point3f v1_cords = v1->P();
				VertexType* mutual_v1 = hedge.V();
				vcg::Point3f mutual_v1_cords = mutual_v1->P();
				hedge.FlipV();
				VertexType* mutual_v2 = hedge.V();
				vcg::Point3f mutual_v2_cords = mutual_v2->P();
				vcg::Point3f center_sphere = (mutual_v1_cords + mutual_v2_cords) / 2;
				double radius = vcg::Distance(center_sphere, mutual_v1_cords);

				hedge.FlipF();
				hedge.FlipE();
				hedge.FlipV();
				VertexType* v2 = hedge.V(); //opposing vertex
				vcg::Point3f v2_cords = v2->P();
				FaceType* opposing_face = hedge.F();
				//print_vertices(current_face);
				//print_vertices(opposing_face);

				// we calculate the distances of the current edge and the relaxed edge
				// we relax the edge if the other edge is shorter than the current one
				double non_mut_dist1 = vcg::Distance(center_sphere, v1_cords);
				double non_mut_dist2 = vcg::Distance(center_sphere, v2_cords);

				if ((non_mut_dist1 <= radius) && (non_mut_dist2 <= radius)) {
					/*
					To relax an edge means, for the two triangles adjacent to
					the edge, to check that each of the two non-mutual vertices
					of these triangles lies outside of the circum-sphere of the
					opposing triangle. If this test fails, the edge is swapped.

					In this case, we know that the distance from the center of the sphere
					to the two non mutual vertices is smaller than the radius, which means that the two
					non mutual vertices are inside the sphere.
					And therefore, the new edge is smaller than the current one.
					*/

					//we change the coordinates of the vertices of each of the triangles
					int i = 0;
					for (i = 0; i < current_face->VN(); i++) {
						if (current_face->V(i) == mutual_v1) {
							current_face->V(i) = v2;
							break;
						}
					}
					assert(i != current_face->VN());
					
					for (i = 0; i < opposing_face->VN(); i++) {
						if (opposing_face->V(i) == mutual_v2) {
							opposing_face->V(i) = v1;
							break;
						}
					}
					assert(i != opposing_face->VN());
					//opposing_face->SetV();
					vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
					//if (!vcg::face::CheckFlipEdge(*current_face, edge_to_swap)) {
					//	std::cout << "edge cant be swapped\n";
					//}
					//vcg::face::FlipEdge(*current_face, edge_to_swap);
					relaxed_edge = true;
					return true;
				}
				else {
					//We don't relax the edge
					return false;
				}
			};

			for (auto fi = mesh.face.begin(); fi != mesh.face.end(); fi++) if(!fi->IsD()) {
				vcg::face::Pos<FaceType> hedge0(&*fi, 0);
				relax_edge(hedge0);

				vcg::face::Pos<FaceType> hedge1(&*fi, 1);
				relax_edge(hedge1);

				vcg::face::Pos<FaceType> hedge2(&*fi, 2);
				relax_edge(hedge2);
			}
			//returns true if any interior edge was relaxed, false otherwise
			vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
			return relaxed_edge;
		};

		while (true) {
			bool s2 = step2();
			//vcg::tri::io::ExporterOFF<MeshType>::Save(mesh, ("step2_" + std::to_string(c) + ".off").c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);

				   //step 3
			if (!s2) break; //no triangles were added, so the patching mesh is complete

				   //step 5
			while (step4()); //while interior edges are relaxed, keep doing step 4
			//otherwise we go back to step 2 (begining of the loop)
			//vcg::tri::io::ExporterOFF<MeshType>::Save(mesh, ("step4_" + std::to_string(c) + ".off").c_str(), vcg::tri::io::Mask::IOM_FACECOLOR);
		}

		//incorrect solutions to fix the color of the faces
		//
		//vcg::tri::UpdateColor<MeshType>::PerFaceConstant(mesh);
		//
		//vcg::tri::UpdateNormal<MeshType>::PerVertex(mesh);
		//vcg::tri::UpdateNormal<MeshType>::PerFace(mesh);
	}
};

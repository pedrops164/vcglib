#pragma once
#include "my_mesh.h"
#include "alg_util.h"
#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
using namespace std;
using namespace vcg;

template<class MeshType>
class RefinementAlg {
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
private:
	double alpha;
	MeshType& m;
	MeshType* new_mesh = nullptr;
	MeshType* aux_mesh = nullptr;
	vector<VertexType*> bV; //border vertices
	vector<FaceType*> bF; //border faces
	vector<double> scale; //vector of scale for each vertex. scale[i] is the scale of the i_th vertexs

public:

	RefinementAlg(MeshType& m, vector<VertexType*> bV, vector<FaceType*> bF): m(m), bV(bV), bF(bF) {
		alpha = sqrt(2);

		aux_mesh = new MeshType();
		tri::RequirePerVertexNormal(*aux_mesh);

		size_t n = bV.size();
		auto ptr_m2 = tri::Allocator<MeshType>::AddVertices(*aux_mesh, n);

		for (auto vi = m.vert.begin(); vi != m.vert.end(); vi++, ptr_m2++) {
			ptr_m2->P() = vi->P();
			ptr_m2->N() = vi->N();
		}
		for (auto fi = m.face.begin(); fi != m.face.end(); fi++) {
			auto v0 = fi->V(0);
			auto v1 = fi->V(1);
			auto v2 = fi->V(2);
			int v0_index = AlgUtil<MeshType>::index_of_vertex(m, v0);
			int v1_index = AlgUtil<MeshType>::index_of_vertex(m, v1);
			int v2_index = AlgUtil<MeshType>::index_of_vertex(m, v2);
			AlgUtil<MeshType>::addTriangle(*aux_mesh,
						AlgUtil<MeshType>::get_vertex(*aux_mesh, v0_index),
						AlgUtil<MeshType>::get_vertex(*aux_mesh, v1_index),
						AlgUtil<MeshType>::get_vertex(*aux_mesh, v2_index));
		}
		tri::UpdateTopology<MeshType>::FaceFace(*aux_mesh);
	}

	MeshType& algorithm() {

			   //step 1
		/*  For each vertex v_i on the hole boundary, compute the
			scale attribute sigma(v_i) as the average length of the
			edges that are adjacent to v_i in the surrounding mesh.
			Initialize the patching mesh as the given hole
			triangulation.
		*/
		for (int i = 0; i < bV.size(); i++) {
			//iterate over all the faces on the border
			//face::Pos<FaceType> cur(bF[i], 0);
			face::Pos<FaceType> cur(bF[i], bV[i]);
			if (!cur.IsBorder()) {
				//if it's not a border, switch the edge (to go to the border)
				cur.FlipE();
				assert(cur.IsBorder());
			}
			Point3f orig = cur.V()->P();
			double sum_distance = 0;
			int count = 0;
			do {
				Point3f dist = cur.VFlip()->P();
				sum_distance += distance(orig, dist);
				count++;
				cur.FlipE();
				cur.FlipF();

			} while (!cur.IsBorder());
			//calculate last edge length
			Point3f dest = cur.VFlip()->P();
			sum_distance += distance(orig, dest);
			count++;

			double avg_distance = sum_distance / count;
			scale.push_back(avg_distance);
			//scale_att.push_back(avg_distance);
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
		//auto step2 = [&]() {
		//	bool triangles_added = false;
		//
		//	vector<FaceType**> faces_pu;
		//	//vector<VertexType**> vertices_pu;
		//	for (auto fi = aux_mesh->face.begin(); fi != aux_mesh->face.end(); fi++) {
		//		VertexType* v0 = fi->V(0);
		//		VertexType* v1 = fi->V(1);
		//		VertexType* v2 = fi->V(2);
		//		Point3f v0_p = v0->P();
		//		Point3f v1_p = v1->P();
		//		Point3f v2_p = v2->P();
		//		int v0_index = index_of_vertex(*aux_mesh, v0);
		//		int v1_index = index_of_vertex(*aux_mesh, v1);
		//		int v2_index = index_of_vertex(*aux_mesh, v2);
		//
		//		Point3f centroid_cords = (v0_p + v1_p + v2_p) / 3;
		//		double scale_v0 = scale[v0_index];
		//		double scale_v1 = scale[v1_index];
		//		double scale_v2 = scale[v2_index];
		//		double scale_centroid = (scale_v0 + scale_v1 + scale_v2) / 3;
		//
		//		//definition of alpha condition
		//		auto check = [&](VertexType* v) {
		//			Point3f dest = v->P();
		//			double dist = distance(centroid_cords, dest);
		//			double value = alpha * dist;
		//			if (value > scale_centroid && value > scale[index_of_vertex(*aux_mesh, v)]) return true;
		//			return false;
		//		};
		//
		//
		//		//For m=i,j,k, if the alpha condition is met, subdivide triangle
		//		if (check(v0) || check(v1) || check(v2)) {
		//			//we add the centroid vertex to the new mesh
		//			// we can only add it once.
		//			triangles_added = true;
		//			//VertexType* centroid = &*tri::Allocator<MeshType>::AddVertices(*aux_mesh, 1);
		//			FaceType* f = &*fi;
		//			faces_pu.push_back(&f);
		//			//centroid->P() = centroid_cords;
		//			scale.push_back(scale_centroid);
		//
		//			//v0 = &aux_mesh->vert[v0_index];
		//			//v1 = &aux_mesh->vert[v1_index];
		//			//v2 = &aux_mesh->vert[v2_index];
		//			//
		//			//fi->V(2) = centroid; //modify current triangle
		//			//
		//			////and add other two triangles
		//			//addTriangle(*aux_mesh, centroid, v0, v2);
		//			//addTriangle(*aux_mesh, centroid, v1, v2);
		//		}
		//		else {
		//			//otherwise we don't add any triangle to the current mesh
		//		}
		//	}
		//
		//	//for each face to be split, we add two faces to the mesh
		//	// because the each face is modified to be one of the three new triangles, and the other two
		//	// are added.
		//	auto f_p = tri::Allocator<MeshType>::AddFaces(*aux_mesh, faces_pu.size() * 2, faces_pu);
		//
		//	//for each face to be added, one vertex is added (the centroid of each face)
		//	auto v_p = tri::Allocator<MeshType>::AddVertices(*aux_mesh, faces_pu.size());
		//
		//	for (FaceType** f : faces_pu) {
		//		FaceType* fi = *f;
		//		VertexType* v0 = fi->V(0);
		//		VertexType* v1 = fi->V(1);
		//		VertexType* v2 = fi->V(2);
		//		Point3f v0_p = v0->P();
		//		Point3f v1_p = v1->P();
		//		Point3f v2_p = v2->P();
		//		int v0_index = index_of_vertex(*aux_mesh, v0);
		//		int v1_index = index_of_vertex(*aux_mesh, v1);
		//		int v2_index = index_of_vertex(*aux_mesh, v2);
		//
		//		Point3f centroid_cords = (v0_p + v1_p + v2_p) / 3;
		//		VertexType* centroid = &*v_p;
		//		v_p++;
		//		centroid->P() = centroid_cords;
		//
		//		fi->V(2) = centroid;
		//		FaceType* new_face1 = &*f_p;
		//		f_p++;
		//		FaceType* new_face2 = &*f_p;
		//		f_p++;
		//		new_face1->V(0) = centroid;
		//		new_face1->V(1) = v0;
		//		new_face1->V(2) = v2;
		//
		//		new_face2->V(0) = centroid;
		//		new_face2->V(1) = v1;
		//		new_face2->V(2) = v2;
		//	}
		//
		//	//after all vertices have been added to new mesh, we update the topology
		//	tri::UpdateTopology<MeshType>::FaceFace(*aux_mesh);
		//	//new_mesh = aux_mesh;
		//	return triangles_added;
		//};

		auto step2 = [&]() {
			bool triangles_added = false;

				   //here create the auxiliar mesh, and copy all vertices of m to it
				   //at the end of this step, we free m, and set it to the new mesh
			new_mesh = new MeshType();
			int n = aux_mesh->VN();
			auto ptr = tri::Allocator<MeshType>::AddVertices(*new_mesh, n);
			for (int i = 0; i < n; i++, ptr++) {
				//update coords and normal of each vertex
				ptr->P() = (&aux_mesh->vert.front())[i].P();
				ptr->N() = (&aux_mesh->vert.front())[i].N();
			}
			for (auto fi = aux_mesh->face.begin(); fi != aux_mesh->face.end(); fi++) {
				auto v0 = fi->V(0);
				auto v1 = fi->V(1);
				auto v2 = fi->V(2);

				int v0_index = AlgUtil<MeshType>::index_of_vertex(*aux_mesh, v0);
				int v1_index = AlgUtil<MeshType>::index_of_vertex(*aux_mesh, v1);
				int v2_index = AlgUtil<MeshType>::index_of_vertex(*aux_mesh, v2);

				Point3f centroid_cords = (v0->P() + v1->P() + v2->P()) / 3;
				double scale_v0 = scale[v0_index];
				double scale_v1 = scale[v1_index];
				double scale_v2 = scale[v2_index];
				double scale_centroid = (scale_v0 + scale_v1 + scale_v2) / 3;

					   //definition of alpha condition
				auto check = [&](VertexType* v) {
					Point3f dest = v->P();
					double dist = distance(centroid_cords, dest);
					double value = alpha * dist;
					if (value > scale_centroid && value > scale[AlgUtil<MeshType>::index_of_vertex(*aux_mesh, v)]) return true;
					return false;
				};

				/*
				 * This lambda function receives a half edge on the current face, and calculates if the edge
				 * needs to be relaxed, and adds the new triangle accordingly.
				 */
				auto replace_relax = [&](face::Pos<FaceType> hedge, VertexType* centroid) {
					face::Pos<FaceType> aux2 = hedge;
					VertexType* mutual_vert1 = aux2.V();
					Point3f mutual_vert1_cords = mutual_vert1->P();
					aux2.FlipV();
					VertexType* mutual_vert2 = aux2.V();
					Point3f mutual_vert2_cords = mutual_vert2->P();
					double current_edge_dist = distance(mutual_vert1_cords, mutual_vert2_cords);
					int mutual_vert1_index = AlgUtil<MeshType>::index_of_vertex(*aux_mesh, mutual_vert1);
					int mutual_vert2_index = AlgUtil<MeshType>::index_of_vertex(*aux_mesh, mutual_vert2);
					//int centroid_index = index_of_vertex(*new_mesh, centroid);
					AlgUtil<MeshType>::addTriangle(*new_mesh, centroid, AlgUtil<MeshType>::get_vertex(*new_mesh, mutual_vert1_index), AlgUtil<MeshType>::get_vertex(*new_mesh, mutual_vert2_index));
					return;
					//if (hedge.IsBorder()) {
					//	//if the half edge is on the border, this means that the adjacent triangle is
					//	// on the surrounding mesh, and we don't want to modify it, so we just add the triangle
					//	addTriangle(*new_mesh, centroid, get_vertex(*new_mesh, mutual_vert1_index), get_vertex(*new_mesh, mutual_vert2_index));
					//	//addTriangle(*new_mesh, centroid, &*(new_mesh->vert.begin()+vert1_index), old2new[mutual_vert2]);
					//	return;
					//}
					//
					//face::Pos<FaceType> aux1 = hedge;
					//aux1.FlipF();
					//aux1.FlipE();
					//aux1.FlipV();
					//FaceType* adj_face = aux1.F();
					////right now, aux is on the non-mutual vertex of the adjacent face
					//VertexType* non_mutual_vert = aux1.V();
					//int non_mutual_vert_index = index_of_vertex(*aux_mesh, non_mutual_vert);
					//Point3f non_mutual_vert_cords = aux1.V()->P();
					//double relaxed_edge_dist = distance(centroid_cords, non_mutual_vert_cords);
					//
					//if (relaxed_edge_dist < current_edge_dist) {
					//	//We relax the edge, because the relaxed edge is smaller than the current one.
					//	//So we add both triangles to the new mesh, and we set the other triangle as visited.
					//	adj_face->SetV();
					//	addTriangle(*new_mesh, &*centroid, get_vertex(*new_mesh, mutual_vert1_index), get_vertex(*new_mesh, non_mutual_vert_index));
					//	addTriangle(*new_mesh, &*centroid, get_vertex(*new_mesh, mutual_vert2_index), get_vertex(*new_mesh, non_mutual_vert_index));
					//
					//}
					//else {
					//	//We don't relax the edge, so we just add the triangle
					//	// the face adjacent to this new triangle remains unchanged, and is yet to be visited.
					//	addTriangle(*new_mesh, centroid, get_vertex(*new_mesh, mutual_vert1_index), get_vertex(*new_mesh, mutual_vert2_index));
					//}
				};

					   //For m=i,j,k, if the alpha condition is met, subdivide triangle and relax outer edges
				if (check(v0) || check(v1) || check(v2)) {
					//we add the centroid vertex to the new mesh
					// we can only add it once.
					triangles_added = true;
					auto centroid = tri::Allocator<MeshType>::AddVertices(*new_mesh, 1);
					centroid->P() = centroid_cords;
					scale.push_back(scale_centroid);

						   //we create the half edge positions on each of the outer edges
					face::Pos<FaceType> hedge0(&*fi, 0, v0);
					face::Pos<FaceType> hedge1(&*fi, 1, v1);
					face::Pos<FaceType> hedge2(&*fi, 2, v2);

					replace_relax(hedge0, &*centroid);
					replace_relax(hedge1, &*centroid);
					replace_relax(hedge2, &*centroid);
				}
				else {
					//otherwise we just add this triangle to the new mesh
					AlgUtil<MeshType>::addTriangle(*new_mesh,
								AlgUtil<MeshType>::get_vertex(*new_mesh, v0_index),
								AlgUtil<MeshType>::get_vertex(*new_mesh, v1_index),
								AlgUtil<MeshType>::get_vertex(*new_mesh, v2_index));
				}
			}
			//after all vertices have been added to new mesh, we update the topology
			tri::UpdateTopology<MeshType>::FaceFace(*new_mesh);

				   //now we free auxiliar mesh, set it to null, and swap the new mesh with the aux mesh
				   //aux_mesh->free(); ???
			aux_mesh = new_mesh;
			return triangles_added;
		};

		//step 4
		/*
		Relax all interior edges of the patching mesh.
		*/
		auto step4 = [&]() {
			//first, set all faces as unvisited
			//for (auto fi = m.face.begin(); fi != m.face.end(); fi++) {
			//	fi->ClearV();
			//}


			bool relaxed_edge = false;

			/*
			 * This lambda function receives a half edge position, and tries to relax the edge. It returns true if the edge
			 * was relaxed, false otherwise. If the half edge is a border, it just returns false.
			 */
			auto relax_edge = [&](face::Pos<FaceType> hedge) {
				if (hedge.IsBorder()) {
					return false; //we don't relax the edge, so we return false
				}

					   //Here we calculate the 4 vertices involved and their coordinates
				auto aux1 = hedge;
				aux1.FlipE();
				aux1.FlipV();
				FaceType* current_face = aux1.F();
				VertexType* v1 = aux1.V(); //internal vertex
				Point3f v1_cords = v1->P();
				VertexType* mutual_v1 = hedge.V();
				Point3f mutual_v1_cords = mutual_v1->P();
				hedge.FlipV();
				VertexType* mutual_v2 = hedge.V();
				Point3f mutual_v2_cords = mutual_v2->P();
				Point3f center_sphere = (mutual_v1_cords + mutual_v2_cords) / 2;
				double radius = distance(center_sphere, mutual_v1_cords);

				hedge.FlipF();
				hedge.FlipE();
				hedge.FlipV();
				VertexType* v2 = hedge.V(); //opposing vertex
				Point3f v2_cords = v2->P();
				FaceType* opposing_face = hedge.F();
				//print_vertices(current_face);
				//print_vertices(opposing_face);

					   // we calculate the distances of the current edge and the relaxed edge
					   // we relax the edge if the other edge is shorter than the current one
				double non_mut_dist1 = distance(center_sphere, v1_cords);
				double non_mut_dist2 = distance(center_sphere, v2_cords);

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
					for (int i = 0; i < current_face->VN(); i++) {
						if (current_face->V(i) == mutual_v1) {
							current_face->V(i) = v2;
							break;
						}
					}
					for (int i = 0; i < opposing_face->VN(); i++) {
						if (opposing_face->V(i) == mutual_v2) {
							opposing_face->V(i) = v1;
							break;
						}
					}
					//opposing_face->SetV();
					relaxed_edge = true;
					tri::UpdateTopology<MeshType>::FaceFace(*new_mesh);
					return true;
				}
				else {
					//We don't relax the edge
					return false;
				}
			};

			for (auto fi = (*new_mesh).face.begin(); fi != (*new_mesh).face.end(); fi++) {
				face::Pos<FaceType> hedge0(&*fi, 0);
				relax_edge(hedge0);

				face::Pos<FaceType> hedge1(&*fi, 1);
				relax_edge(hedge1);

				face::Pos<FaceType> hedge2(&*fi, 2);
				relax_edge(hedge2);
			}
			//returns true if any interior edge was relaxed, false otherwise
			tri::UpdateTopology<MeshType>::FaceFace(*new_mesh);
			return relaxed_edge;
		};


		while (true) {
			cout << aux_mesh->FN() << " " << aux_mesh->VN() << endl;
			bool s2 = step2();
			cout << new_mesh->FN() << " " << new_mesh->VN() << endl;

				   //step 3
			if (!s2) break; //no triangles were added, so the patching mesh is complete

				   //step 5
			while (step4()) cout << "cp0" << endl; //while interior edges are relaxed, keep doing step 4
			//otherwise we go back to step 2 (begining of the loop)
		}
		return *new_mesh;
	}
};

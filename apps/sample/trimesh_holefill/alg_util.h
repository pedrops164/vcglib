#pragma once

#include "my_mesh.h"

template<class MeshType>
class AlgUtil {
public:
	typedef typename MeshType::VertexType VertexType;
	typedef typename MeshType::FaceType FaceType;

	/*
	 * Returns the index of the vertex
	 * Assumes the vertex is in the mesh
	 * Returns -1 if the vertex doesnt belong to the mesh
	 */
	static int index_of_vertex(MeshType& mesh, VertexType* v) {
		assert(mesh.vert.size() > 0); //is this needed?
		auto index_front = &mesh.vert.front();
		auto index_back = &mesh.vert.back();
		if (v < index_front || v > index_back) {
			return -1;
		}
		//assert(v >= index_front && v <= index_back);
		return v - index_front;
	}

	static VertexType* get_vertex(MeshType& mesh, int index) {
		assert(mesh.vert.size() > index);
		return (VertexType*) &(mesh.vert.begin()[index]); //careful here
	}

	//static void print_vertices(MeshType& mesh) {
	//	for (auto fi = mesh.vert.begin(); fi != mesh.vert.end(); fi++) {
	//		cout << fi->P().X() << ", " << fi->P().Y() << ", " << fi->P().Z() << endl;
	//	}
	//}
	//static void print_vertices(FaceType* fi) {
	//	cout << "printing vertices" << endl;
	//	for (int i = 0; i < fi->VN(); i++) {
	//		cout << fi->V(i)->P().X() << ", " << fi->V(i)->P().Y() << ", " << fi->V(i)->P().Z() << endl;
	//	}
	//}

	/*
	 * Receives a reference to a mesh, and returns a mesh with all vertices and faces of the receiving mesh
	 */
	static MeshType* mesh_copy(MeshType& m) {
		MeshType* ret = new MeshType();
		int n_vert = m.VN();
		auto ptr_verts = vcg::tri::Allocator<MeshType>::AddVertices(*ret, n_vert);
		//adds all vertices
		for (int i = 0; i < n_vert; i++, ptr_verts++) {
			//update coords and normal of each vertex
			ptr_verts->P() = m.vert[i].P();
			ptr_verts->N() = m.vert[i].N();
		}

			   //adds all faces
		for (auto fi = m.face.begin(); fi != m.face.end(); fi++) {
			auto v0 = fi->V(0);
			auto v1 = fi->V(1);
			auto v2 = fi->V(2);
			int v0_index = index_of_vertex(m, v0);
			int v1_index = index_of_vertex(m, v1);
			int v2_index = index_of_vertex(m, v2);
			vcg::tri::Allocator<MeshType>::AddFace(*ret,
						get_vertex(*ret, v0_index),
						get_vertex(*ret, v1_index),
						get_vertex(*ret, v2_index));
		}
		return ret;
	}

	/*
	 * Frees a mesh
	 */
	//static void free_mesh(MeshType& m) {
	//	for (auto fi = m.face.begin(); fi != m.face.end(); fi++) {
	//		tri::Allocator<MeshType>::DeleteFace(m, *fi);
	//	}
	//	for (auto vi = m.vert.begin(); vi != m.vert.end(); vi++) {
	//		tri::Allocator<MeshType>::DeleteVertex(m, *vi);
	//	}
	//}

};


#pragma once

#include <unordered_map>
#include <vector>
#include <vcg/complex/algorithms/hole.h>

template<class MeshType>
class TriangulationAlg {
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
	typedef typename vcg::face::Pos<FaceType> PosType;
	typedef typename vcg::tri::Hole<MeshType>::Weight Weight;
private:
	MeshType& mesh;
	std::vector<std::vector<Weight>> w_matrix;
	std::vector<std::vector<int>> delta_matrix;
	std::vector<PosType> vv;

	/*
	 * The trace function recovers the triangulation by accessing the matrix
	 * It is called recursively
	 */
	void trace(int i, int k) {
		auto front = mesh.vert.begin();
		if (i + 2 == k) {
			// In this case, we just append the triangle to the mesh.
			vcg::tri::Allocator<MeshType>::AddFace(mesh, &front[i], &front[i + 1], &front[i + 2]);
		}
		else {
			int o = delta_matrix[i][k];
			if (o != i + 1) trace(i, o);
			vcg::tri::Allocator<MeshType>::AddFace(mesh, &front[i], &front[o], &front[k]);
			if (o != k - 1) trace(o, k);
		}
	}
public:

	TriangulationAlg(MeshType& mesh, std::vector<PosType> vv): mesh(mesh), vv(vv) {
		vcg::tri::RequirePerVertexNormal(mesh);
		size_t n = vv.size();
		auto ptr = vcg::tri::Allocator<MeshType>::AddVertices(mesh, n);
		for (int i = 0; i < n; i++, ptr++) {
			//update coords and normal of each vertex
			ptr->P() = vv[i].v->P();
			ptr->N() = vv[i].v->N();
		}
	}

	//runs the triangulation algorithm
	void algorithm() {
		//hole size
		int n = vv.size();

		w_matrix.clear();
		w_matrix.resize(n, std::vector<Weight>(n, Weight()));

		delta_matrix.resize(n, std::vector<int>(n, 0));

		//step 1

		for (int i = 0; i < n - 1; i++) {
			w_matrix[i][i + 1] = Weight(0, 0);
		}

		for (int i = 0; i < n - 2; i++) {
			//calculate the angle seperately
			//in step 1, the angle(i,i+1,i+2) is calculated against the triangles adjacent to edges
			// (i,i+1) and (i+1, i+2) in the surrounding mesh.

			w_matrix[i][i+2] = vcg::tri::Hole<MeshType>::computeWeight(i, i + 1, i + 2, vv, delta_matrix);
		}

		int j = 2;
		do {
			//step 2
			j++;
			for (int i = 0; i < n - j; i++) {
				int k = i + j;
				Weight min_value;
				int min_m = -1;
				for (int m = i + 1; m < k; m++) {
					Weight a = w_matrix[i][m];
					Weight b = w_matrix[m][k];
					Weight newval = a + b + vcg::tri::Hole<MeshType>::computeWeight(i, m, k, vv, delta_matrix);

					if (newval < min_value) {
						min_value = newval;
						min_m = m;
					}
				}
				delta_matrix[i][k] = min_m;
				w_matrix[i][k] = min_value;
			}

		} while (j < n - 1); // step 3 : go back to step 2 or not

		//Step 4
		// Recover the triangulation
		// In this case we are going to add the triangles to the mesh
		trace(0, n - 1);
		vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
	}
};


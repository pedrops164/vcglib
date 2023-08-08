#pragma once
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
#include<Eigen/SparseQR>

template<class MeshType>
class FairingAlg {
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
private:
	int bv_n;
	MeshType& mesh;
	std::vector<int> internal_verts; //vector of internal vertices (vertices not on the border)

	/*
	 * Returns a vector with the neighbors vertices of v
	 * v must belong to the mesh
	 * But the function might returns vertices that don't belong to the mesh (in case v is on the border)
	 * 
	 * Assumes there is VF Adjacency in the mesh
	 *
	 */
	std::vector<VertexType*> get_neighbors(VertexType* v) {
		//int v_idx = AlgUtil<MeshType>::index_of_vertex(mesh, v);
		int v_idx = vcg::tri::Index(mesh, v);
		assert(v_idx >= 0); //makes sure v belongs to mesh
		std::vector<VertexType*> ret;
		vcg::face::VVStarVF<FaceType>(v, ret);
		//if (v_idx < bv_n && v_idx >= 0) {
		if (v->IsB()) {
			//vertex is on the border. We also need to add the vertices on the surrounding mesh.

				   //we establish a halfedge position, on a face and vertex which are on the border
			assert(bV[v_idx] != v); //assert that they are different vertices
			assert(bV[v_idx]->P() == v->P()); //assert that the coordinates are the same
			vcg::face::Pos<FaceType> start(bF[v_idx], bV[v_idx]);

			//the halfedge might be on the border, so we switch the edge in this case
			if (start.IsBorder()) {
				start.FlipE();
			}

			/*
			 * In this cycle, we add all neighbor vertices which are not on the border, because the ones
			 * on the border were already added in vcg::edge::VVStarVF(v, ret);
			 * So we don't want to add duplicate vertices (they are on different meshes, but they are on the same position).
			 */
			do {
				ret.push_back(start.VFlip());
				start.FlipF();
				start.FlipE();
			} while (!start.IsBorder());

		}
		return ret;
	}
	std::vector<VertexType*> bV;
	std::vector<FaceType*> bF;
public:

	FairingAlg(MeshType& mesh, std::vector<VertexType*> bV, std::vector<FaceType*> bF) : mesh(mesh), bV(bV), bF(bF) {
		bv_n = bV.size();
		//adds to internal_verts the indices of the vertices which are not on the border.
		for (int i = 0; i < mesh.VN(); i++) {
			if (i >= bv_n) internal_verts.push_back(i); //if the vertex is not on the border, add to internal_verts
		}
		vcg::tri::UpdateTopology<MeshType>::VertexFace(mesh);
		vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
		vcg::tri::UpdateFlags<MeshType>::FaceBorderFromFF(mesh);
	}

	static float umbrella_operator(VertexType* v1, VertexType* v2) {
		// uniform umbrella-operator
		return 1;

		// scale-dependent umbrella-operator
		//return vcg::Distance(v1->P(), v2->P());

		//harmonic umbrella-operator
		// ......
	}

	/*
	 * Gets a half edge on the edge whose umbrella operator needs to be calculated
	 * Assumes the hedge has 2 adjacent triangles!!!
	 */
	//float harmonic(vcg::face::Pos<FaceType> hedge) {
	//	assert(hedge.F() != hedge.FFlip()); //ensures the hedge has 2 adjacent triangles
	//	vcg::face::Pos<FaceType> vk1 = hedge;
	//	vcg::face::Pos<FaceType> vk2 = hedge;
	//	VertexType* vi = hedge.V();
	//	VertexType* vj = hedge.VFlip();
	//	vk1.FlipE();
	//	vk1.FlipV();
	//	float angle_vk1 = vk1.AngleRad(); //angle in radians
	//	vk2.FlipF();
	//	vk2.FlipE();
	//	vk2.FlipV();
	//	float angle_vk2 = vk2.AngleRad(); //angle in radians
	//
	//		   //cos and sin functions take radians
	//	return (cos(angle_vk1) / sin(angle_vk1)) + (cos(angle_vk2) / sin(angle_vk2));
	//	//cos(a) / sin(a) = cot(a)
	//}

	/*
	 * For each neighbor of v, sums to the result the value of cu(v, v_i)
	 * For the uniform operator, cu(v) just returns the number of neighbors of v
	 */
	float cu(VertexType* v) {
		//assert(!v->IsB());
		float cu_v = 0;
		std::vector<VertexType*> neighbors = get_neighbors(v);
		for (VertexType* v_i : neighbors) {
			cu_v += umbrella_operator(v, v_i);
		}
		return cu_v;
	}

	void algorithm() {

		typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;
		typedef Eigen::VectorXd Vector;

		int n = internal_verts.size();
		int total_n = mesh.vert.size();
		SparseMatrix sparseMatrix(n, n);
		SparseMatrix vect(n, 3);

		for (int i = 0; i < n; i++) {
			for (int col = 0; col < 3; col++) {
				vect.insert(i, col) = 0;
			}
		}

		for (int j = 0; j < n; j++) {
			int index = internal_verts[j];

			std::vector<double> points(n, 0);

			VertexType* v = &mesh.vert[index];
			std::vector<VertexType*> v_i_vec = get_neighbors(v);
			float cu_v = cu(v);
			for (VertexType* v_i : v_i_vec) {
				int v_i_idx = vcg::tri::Index(mesh, v_i);
				//curr_v is the vertex on the other end of the edge
				float value1 = 2 * (1 / -cu_v) * umbrella_operator(v, v_i);
				assert(v_i != v);
				if (v_i_idx >= bv_n && v_i_idx < mesh.VN()) {
					//the vertex v_j is not on the border
					points[v_i_idx - bv_n] += value1;
				}
				else {
					//the values of the vertices are added to the right side of the equation
					vect.coeffRef(j, 0) = vect.coeffRef(j, 0) - value1 * (double)v_i->P().X();
					vect.coeffRef(j, 1) = vect.coeffRef(j, 1) - value1 * (double)v_i->P().Y();
					vect.coeffRef(j, 2) = vect.coeffRef(j, 2) - value1 * (double)v_i->P().Z();
				}
				std::vector<VertexType*> v_j_vec = get_neighbors(v_i);
				float cu_v_i = cu(v_i);
				for (VertexType* v_j : v_j_vec) {
					int v_j_idx = vcg::tri::Index(mesh, v_j);
					float multiplier = umbrella_operator(v_i, v_j) / (cu_v * cu_v_i);
					//if the vertex shouldnt be modified (is on the border, or doesnt belong to the mesh), then
					// the value to be added to the vertex is set to 0 (it isn't added to the left side of the equation)
					//if (v_j != v)
					if (v_j_idx >= bv_n && v_j_idx < mesh.VN()) {
						//the vertex v_j is not on the border
						points[v_j_idx-bv_n] += multiplier;
					}
					else {
						//the values of the vertices are added to the right side of the equation
						vect.coeffRef(j, 0) = vect.coeffRef(j, 0) - multiplier *(double)v_j->P().X();
						vect.coeffRef(j, 1) = vect.coeffRef(j, 1) - multiplier *(double)v_j->P().Y();
						vect.coeffRef(j, 2) = vect.coeffRef(j, 2) - multiplier * (double)v_j->P().Z();
					}
				}
			}
			points[index - bv_n] += 1;

			for (int i = 0; i < n; i++) {
				sparseMatrix.insert(j, i) = points[i];
			}
		}

		//std::cout << "printing linear system! " << std::endl;
		//for (int i = 0; i < n; i++) {
		//	//left side of the equation
		//	for (int j = 0; j < n; j++) {
		//		std::cout << sparseMatrix.coeffRef(i, j) << " ";
		//	}
		//	std::cout << "| ";
		//	//right side of the equation
		//	for (int col = 0; col < 3; col++) {
		//		std::cout << vect.coeffRef(i, col) << " ";
		//	}
		//	std::cout << std::endl;
		//}

		//Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> solver;
		Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver;
		
		sparseMatrix.makeCompressed();
		solver.compute(sparseMatrix);
		if (solver.info() != Eigen::Success) {
			// We calculate the information and check if it is a success or not
			std::cout << "error compute" << std::endl;
			//return mesh;
		}

		SparseMatrix sol = solver.solve(vect); //solution computation
		if (solver.info() != Eigen::Success) {
			// Here the solving has failed
			std::cout << "error solving solution" << std::endl;
			//return mesh;
		}

		//std::cout << "printing solution:\n";
		//for (int i = 0; i < n; i++) {
		//	//left side of the equation
		//	for (int j = 0; j < 3; j++) {
		//		std::cout << sol.coeffRef(i, j) << " ";
		//	}
		//	std::cout << std::endl;
		//}

		for (int i = 0; i < internal_verts.size(); i++) {
			VertexType* v = &mesh.vert[i+bv_n];
			v->P() = vcg::Point3f(sol.coeffRef(i, 0), sol.coeffRef(i, 1), sol.coeffRef(i, 2));
		}
	}

};
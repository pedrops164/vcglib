#pragma once
#include "my_mesh.h"
#include "alg_util.h"
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
#include<Eigen/SparseQR>
#include <unordered_set>

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
	 * Assumes there is VV adjecency on the mesh
	 */
	std::vector<VertexType*> get_neighbors(VertexType* v) {
		int v_idx = AlgUtil<MeshType>::index_of_vertex(mesh, v);
		assert(v_idx >= 0); //makes sure v belongs to mesh
		std::vector<VertexType*> ret;
		vcg::edge::VVStarVE(v, ret);
		if (v_idx < bv_n) {
			//vertex is on the border. We also need to add the vertices on the surrounding mesh.

				   //we establish a halfedge position, on a face and vertex which are on the border
			vcg::face::Pos<FaceType> start(bF[v_idx], bV[v_idx]);

			//the halfedge might be on the border, so we switch the edge in this case
			if (start.IsBorder()) {
				start.FlipE();
			}

			/*
			 * In this cycle, we add all neighbor vertices which are not on the border, because the ones
			 * on the border were already added in vcg::edge::VVStarVE(v, ret);
			 * So we don't want to add duplicate vertices (they are on different meshes, but they are on the same position).
			 */
			do {
				ret.push_back(start.VFlip());
				start.FlipF();
				start.FlipE();
			} while (!start.IsBorder());
			//ret.push_back(start.VFlip());

		}
		if (v_idx == 9) {
			for (auto fi = mesh.face.begin(); fi != mesh.face.end(); ++fi) {
				if (fi->V(0) == v || fi->V(1) == v || fi->V(2) == v) fi->C() = vcg::Color4b::Blue;
			}
			std::cout << ret.size() << " <---------" << std::endl;
			for (VertexType* v2 : ret) {
				std::cout << v2->P().X() << ", " << v2->P().Y() << ", " << v2->P().Z() << "\n";
			}
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
		for (int idx : internal_verts) {
			std::cout << idx << " ";
		}
		std::cout << std::endl;
		vcg::tri::UpdateTopology<MeshType>::AllocateEdge(mesh);
		vcg::tri::UpdateTopology<MeshType>::VertexEdge(mesh);
	}

	//float scale_dependent(VertexType* v1, VertexType* v2) {
	//	return vcg::Distance(v1->P(), v2->P());
	//}

	float uniform(VertexType* v1, VertexType* v2) {
		return 1;
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
		assert(!v->IsB());
		float cu_v = 0;
		std::vector<VertexType*> vert_vec;
		vcg::edge::VVStarVE(v, vert_vec);
		if (vert_vec.size() == 0) std::cout << "is empty" << std::endl;

		for (VertexType* curr_v : vert_vec) {
			assert(curr_v != v);
			cu_v += uniform(v, curr_v);
		}
		return cu_v;
	}

	vcg::Point3f U_cu(VertexType* v) {
		assert(!v->IsB());

			   //cu(v) calculation
		float cu_v = cu(v);

		std::vector<VertexType*> vert_vec;
		vcg::edge::VVStarVE(v, vert_vec);

		vcg::Point3f ret(0, 0, 0);
		for (VertexType* curr_v : vert_vec) {
			//curr_v is the vertex on the other end of the edge
			assert(curr_v != v);
			ret += uniform(v, curr_v) * curr_v->P();
			//cout << "cords - " << curr_v->P().X() << ", " << curr_v->P().X() << ", " << curr_v->P().X() << endl;
		}
		//cout << endl;

		ret *= (1 / cu_v);
		ret += -(v->P());
		return ret;
	}

	vcg::Point3f U2_cu(VertexType* v) {
		assert(!v->IsB());
		float cu_v = cu(v);
		vcg::Point3f u_cu_v = U_cu(v);

		std::vector<VertexType*> vert_vec;
		vcg::edge::VVStarVE(v, vert_vec);

		vcg::Point3f ret(0, 0, 0);
		for (VertexType* curr_v : vert_vec) {
			//curr_v is the vertex on the other end of the edge
			assert(curr_v != v);
			ret += (uniform(v, curr_v) * U_cu(curr_v));
		}
		ret *= (1 / cu_v);

		ret = ret - u_cu_v;

		return ret;
	}

	MeshType& algorithm() {

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
			//float cu_v = cu(v);
			float cu_v = v_i_vec.size();
			for (VertexType* v_i : v_i_vec) {
				//float cu_v_i = cu(v_i);
				//curr_v is the vertex on the other end of the edge
				assert(v_i != v);
				std::vector<VertexType*> v_j_vec = get_neighbors(v_i);
				float cu_v_i = v_j_vec.size();
				for (VertexType* v_j : v_j_vec) {
					int v_j_idx = AlgUtil<MeshType>::index_of_vertex(mesh, v_j);
					float multiplier = (1 / -cu_v) * (1 / cu_v_i);
					//if the vertex shouldnt be modified (is on the border, or doesnt belong to the mesh), then
					// the value to be added to the vertex is set to 0 (it isn't added to the left side of the equation)
					//if (v_j != v)
					if (v_j_idx >= bv_n) {
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

		std::cout << "printing linear system! " << std::endl;
		for (int i = 0; i < n; i++) {
			//left side of the equation
			for (int j = 0; j < n; j++) {
				std::cout << sparseMatrix.coeffRef(i, j) << " ";
			}
			std::cout << "| ";
			//right side of the equation
			for (int col = 0; col < 3; col++) {
				std::cout << vect.coeffRef(i, col) << " ";
			}
			std::cout << std::endl;
		}

		//Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> solver;
		//solver.setTolerance(0);
		//solver.setMaxIterations(100);
		Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver;
		//Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> solver;
		//Eigen::BiCGSTAB<SparseMatrix> solver;
		
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
		//std::cout << "#iterations:     " << solver.iterations() << std::endl;
		//std::cout << "estimated error: " << solver.error() << std::endl;

		std::cout << "printing solution:\n";
		for (int i = 0; i < n; i++) {
			//left side of the equation
			for (int j = 0; j < 3; j++) {
				std::cout << sol.coeffRef(i, j) << " ";
			}
			std::cout << std::endl;
		}

		for (int i = 0; i < internal_verts.size(); i++) {
			VertexType* v = &mesh.vert[i+bv_n];
			v->P() = vcg::Point3f(sol.coeffRef(i, 0), sol.coeffRef(i, 1), sol.coeffRef(i, 2));
		}

		return mesh;
	}

};

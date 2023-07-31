
#pragma once

#include "my_mesh.h"
#include <unordered_map>
#include <vector>
#include "map_EF.h"
#include "alg_util.h"


struct AngleArea {
	AngleArea() : angle(0), area(0) {} //fills with dummy values
	AngleArea(double angle, double area) : angle(angle), area(area) {}

	double angle;
	double area;

	AngleArea operator+(const AngleArea& other) const {
		return AngleArea(std::max(angle, other.angle), area + other.area);
	}

	bool operator<(const AngleArea& other) const
	{
		// first compares angle, and then compares area
		return std::tie(angle, area) < std::tie(other.angle, other.area);
	}
};


template<class MeshType>
class TriangulationAlg {
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
private:
	MeshType mesh;
	EdgeFaceMap<MeshType> eFM;
	std::vector<VertexType*> bV;
	std::vector<std::vector<AngleArea>> w_matrix;
	std::vector<std::vector<int>> delta_matrix;


	/*
	 * Calculates the dihedral angle between two triangles
	 * How to calculate normal of a triangle? Maybe there's a method to calculate directly
	 */
	static double calculateDihedralAngle(
		VertexType* f1_v0,
		VertexType* f1_v1,
		VertexType* f1_v2,
		VertexType* f2_v0,
		VertexType* f2_v1,
		VertexType* f2_v2)
	{

			   //edge vectors of first triangle
		vcg::Point3f f1_vec1 = f1_v1->P() - f1_v0->P();
		vcg::Point3f f1_vec2 = f1_v2->P() - f1_v1->P();

			   //normal of the first triangle
		vcg::Point3f normal1(
			f1_vec1.Y() * f1_vec2.Z() - f1_vec1.Z() * f1_vec2.Y(),
			f1_vec1.Z() * f1_vec2.X() - f1_vec1.X() * f1_vec2.Z(),
			f1_vec1.X() * f1_vec2.Y() - f1_vec1.Y() * f1_vec2.X()
			);

			   //edge vectors of second triangle
		vcg::Point3f f2_vec1 = f2_v1->P() - f2_v0->P();
		vcg::Point3f f2_vec2 = f2_v2->P() - f2_v1->P();

			   //normal of the second triangle
		vcg::Point3f normal2(
			f2_vec1.Y() * f2_vec2.Z() - f2_vec1.Z() * f2_vec2.Y(),
			f2_vec1.Z() * f2_vec2.X() - f2_vec1.X() * f2_vec2.Z(),
			f2_vec1.X() * f2_vec2.Y() - f2_vec1.Y() * f2_vec2.X()
			);

			   // Get the face normals
			   //const vcg::Point3f& normal1 = face1->N(); //dummy normal???
			   //const vcg::Point3f& normal2 = face2->N(); //dummy normal???

			   // Calculate the magnitudes of the face normals
		double magnitude1 = sqrt(normal1.X() * normal1.X() + normal1.Y() * normal1.Y() + normal1.Z() * normal1.Z());
		double magnitude2 = sqrt(normal2.X() * normal2.X() + normal2.Y() * normal2.Y() + normal2.Z() * normal2.Z());

			   // Normalize the face normals
		vcg::Point3f n1(normal1.X() / magnitude1, normal1.Y() / magnitude1, normal1.Z() / magnitude1);
		vcg::Point3f n2(normal2.X() / magnitude2, normal2.Y() / magnitude2, normal2.Z() / magnitude2);

			   // Calculate the dot product between the normalized face normals
		double dotProduct = n1.X() * n2.X() + n1.Y() * n2.Y() + n1.Z() * n2.Z(); //????

			   // Calculate the dihedral angle in radians
		double dihedralAngle = acos(dotProduct);

			   // Convert the angle to degrees if desired
			   // double dihedralAngleDegrees = dihedralAngle * 180.0 / M_PI;

		return dihedralAngle;
	}
	/*
	 * The trace function recovers the triangulation by accessing the matrix
	 * It is called recursively
	 */
	void trace(int i, int k) {
		auto front = mesh.vert.begin();
		if (i + 2 == k) {
			// In this case, we just append the triangle to the mesh.
			vcg::tri::Allocator<MeshType>::AddFace(mesh, &front[i], &front[i+1], &front[i+2]);
		}
		else {
			int o = getDelta(i, k);
			if (o != i + 1) trace(i, o);
			vcg::tri::Allocator<MeshType>::AddFace(mesh, &front[i], &front[o], &front[k]);
			if (o != k - 1) trace(o, k);
		}
	}
public:

	TriangulationAlg(std::vector<VertexType*> bV, EdgeFaceMap<MeshType> eFM): bV(bV), eFM(eFM) {
		vcg::tri::RequirePerVertexNormal(mesh);
		size_t n = bV.size();
		auto ptr = vcg::tri::Allocator<MeshType>::AddVertices(mesh, n);
		for (int i = 0; i < n; i++, ptr++) {
			//update coords and normal of each vertex
			ptr->P() = bV[i]->P();
			ptr->N() = bV[i]->N();
		}

		//initializes the matrix as n*n, filled with dummy values.
		w_matrix = std::vector<std::vector<AngleArea>>(n, std::vector<AngleArea>(n, AngleArea()));
		delta_matrix = std::vector<std::vector<int>>(n, std::vector<int>(n, 0));
	}

	int getDelta(int i, int k) {
		return delta_matrix[k][i];
	}

	// sets the delta value of (i,k)
	// stores in the bottom part of the matrix, so that it isn't overwritten with the W values.
	void setDelta(int i, int k, int value) {
		delta_matrix[k][i] = (value);
	}

	//defines the value of W_i,j
	void setW(int i, int j, double angle, double area) {
		w_matrix[i][j] = AngleArea(angle, area);
	}

	void setW(int i, int j, AngleArea anglearea) {
		w_matrix[i][j] = anglearea;
	}

	AngleArea getW(int i, int j) {
		return w_matrix[i][j];
	}


	//runs the triangulation algorithm (first version)
	MeshType& algorithm() {
		int border_vertex_count = bV.size();

			   //step 1

		for (int i = 0; i < border_vertex_count - 1; i++) {
			setW(i, i + 1, 0, 0);
		}

		for (int i = 0; i < border_vertex_count - 2; i++) {
			//calculate the angle seperately
			//in step 1, the angle(i,i+1,i+2) is calculated against the triangles adjacent to edges
			// (i,i+1) and (i+1, i+2) in the surrounding mesh.
			// we can create a position in the vertex i+1, and spin around it, looking for triangles
			// that are next to a border, appending these triangles to a list,
			// and then calculating the max dihedral angle.

				   //cout << bV[i] << "," << bV[i + 1] << "," << bV[i + 2] << endl;
			FaceType* adj_face_1 = eFM.get(bV[i], bV[i + 1]);
			FaceType* adj_face_2 = eFM.get(bV[i + 1], bV[i + 2]);
			setW(i, i + 2, omega(bV[i], bV[i+1], bV[i+2],
								 { {adj_face_1->V(0), adj_face_1->V(1), adj_face_1->V(2)},
								  {adj_face_2->V(0), adj_face_2->V(1), adj_face_2->V(2)}}));
		}

		int j = 2;
		int n = bV.size();
		do {
			FaceType* surr_face = eFM.get(bV[0], bV[n - 1]);
			//step 2
			j++;
			for (int i = 0; i < n - j; i++) {
				int k = i + j;
				int m = i + 1;
				std::vector<std::vector<VertexType*>> adj_faces =
					{ {bV[i], bV[getDelta(i,m)], bV[m]},
					 {bV[m], bV[getDelta(m,k)], bV[k]} };
				if (i == 0 && k == n - 1) {
					adj_faces.push_back({ surr_face->V(0), surr_face->V(1), surr_face->V(2) });
				}
				AngleArea min_value = getW(i, m) + getW(m, k) + omega(bV[i], bV[m], bV[k], adj_faces);
				int min_m = m;
				for (m = m + 1; m < k; m++) {
					adj_faces = { {bV[i], bV[getDelta(i,m)], bV[m]},
								 {bV[m], bV[getDelta(m,k)], bV[k]} };
					if (i == 0 && k == n - 1) {
						FaceType* surr_face = eFM.get(bV[0], bV[n - 1]);
						adj_faces.push_back({ surr_face->V(0), surr_face->V(1), surr_face->V(2) });
					}
					AngleArea current = getW(i, m) + getW(m, k) + omega(bV[i], bV[m], bV[k], adj_faces);
					if (current < min_value) {
						min_value = current;
						min_m = m;
					}
				}
				setW(i, k, min_value);
				setDelta(i, k, min_m);
			}

		} while (j < n - 1); // step 3 : go back to step 2 or not

			   //Step 4
			   // Recover the triangulation
			   // In this case we are going to add the triangles to the mesh
		trace(0, n - 1);
		vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);

		return mesh;
	}


	/*
	 * Calculates the omega function of the triangulation algorithm
	 * For now, just calculates area between the three vertices
	 *
	 * Now it calculates the AngleArea struct of the current triangle, taking into account the
	 * maximum dihedral angle between the current triangle and the adjacent triangles.
	 */
	static AngleArea omega(VertexType* v1, VertexType* v2, VertexType* v3, std::vector<std::vector<VertexType*>> adj_faces) {
		double area = omega_area(v1, v2, v3);
		assert(adj_faces.size() > 0);
		double max_dihedral_angle = 0;
		for (std::vector<VertexType*> f : adj_faces) {
			VertexType* other_v1 = f[0];
			VertexType* other_v2 = f[1];
			VertexType* other_v3 = f[2];
			double angle = calculateDihedralAngle(v1, v2, v3, other_v1, other_v2, other_v3);
			if (angle > max_dihedral_angle) max_dihedral_angle = angle;
		}
		return AngleArea(max_dihedral_angle, area);
	}

	static double omega_area(VertexType* v1, VertexType* v2, VertexType* v3) {
		vcg::Point3f edge1 = v2->P() - v1->P();
		vcg::Point3f edge2 = v3->P() - v1->P();

			   // Calculate the cross product
		vcg::Point3f crossProduct;
		crossProduct[0] = edge1[1] * edge2[2] - edge1[2] * edge2[1];
		crossProduct[1] = edge1[2] * edge2[0] - edge1[0] * edge2[2];
		crossProduct[2] = edge1[0] * edge2[1] - edge1[1] * edge2[0];

			   // Calculate the area as half the length of the cross product
		double area = 0.5 * crossProduct.Norm();

		return area;
	}

};


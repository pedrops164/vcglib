#pragma once

#include <unordered_map>
#include <vector>
#include "my_mesh.h"
using namespace std;

template<class MeshType>
struct EdgeKey
{
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
	EdgeKey(VertexType* first, VertexType* second) : first(first), second(second) {}

	VertexType* first;
	VertexType* second;

	bool operator==(const EdgeKey& other) const
	{
		return (first == other.first
			&& second == other.second) ||
			(first == other.second
				&& second == other.first);
	}
};

//template <>
template<class MeshType>
struct std::hash<EdgeKey<MeshType>>
{
	std::size_t operator()(const EdgeKey<MeshType>& k) const
	{
		using std::size_t;
		using std::hash;
		using std::string;

		// Compute individual hash values for first,
		// second and third and combine them using XOR
		// and bit shifting:

		return (size_t)k.first + (size_t)k.second; //not good hash function?
	}
};

template<class MeshType>
class EdgeFaceMap {
public:
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
	//could instead map edges to faces. This way from an edge on the border I can get its face
	//this way, it maps starting vertex to end vertex of the edge, to the face on the border.
	unordered_map<EdgeKey<MeshType>, FaceType*> edgeFaceMap;

	/*
	* Adds a face to the vector associated with the vertex pointer. Makes sure that the face
	* hasn't been added already!
	*/
	void add(VertexType* start, VertexType* end, FaceType* face) {
		edgeFaceMap[EdgeKey<MeshType>(start, end)] = face;
	}

	/*
	 * Returns the face associated with the given edge on the border.
	 */
	FaceType* get(VertexType* start, VertexType* end){

		auto face = edgeFaceMap.find(EdgeKey<MeshType>(start, end));
		if (face == edgeFaceMap.end()) {
			cout << "couldnt find the face of this edge" << endl;
			return nullptr;
		}
		return face->second;
	}
};


#pragma once
#include<vcg/complex/complex.h>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
    vcg::Use<MyEdge>     ::AsEdgeType,
    vcg::Use<MyFace>     ::AsFaceType> {};

class MyVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags, vcg::vertex::VEAdj> {};
class MyFace : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj, vcg::face::VertexRef, vcg::face::Color4b, vcg::face::BitFlags > {};
class MyEdge : public vcg::Edge<   MyUsedTypes, vcg::edge::VEAdj, vcg::edge::VertexRef> {};

class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace>, std::vector<MyEdge>  > {};
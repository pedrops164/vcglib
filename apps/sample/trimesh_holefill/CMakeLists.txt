cmake_minimum_required(VERSION 3.13)
project(trimesh_holefill)

if (VCG_HEADER_ONLY)
	set(SOURCES
		alg_util.cpp
		alg_util.h
		fairing_alg.cpp
		fairing_alg.h
		map_EF.cpp
		map_EF.h
		my_mesh.h
		refinement_alg.cpp
		refinement_alg.h
		triangulation_alg.cpp
		triangulation_alg.h
		trimesh_holefill.cpp
		trimesh_holefill.h
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_holefill
	${SOURCES})

target_link_libraries(
	trimesh_holefill
	PUBLIC
		vcglib
	)

cmake_minimum_required(VERSION 3.13)
project(lagrange_cpp)

set(CMAKE_CXX_STANDARD 14)

link_directories(/usr/local/lib)
include_directories(/usr/local/include)
include_directories(src)

add_executable(lagrange_cpp
        src/AncSplit.cpp
        src/AncSplit.h
        src/BayesianBioGeo.cpp
        src/BayesianBioGeo.h
        src/BayesianBioGeoAllDispersal.cpp
        src/BayesianBioGeoAllDispersal.h
        src/BioGeoTree.cpp
        src/BioGeoTree.h
        src/BioGeoTreeTools.cpp
        src/BioGeoTreeTools.h
        src/BranchSegment.cpp
        src/BranchSegment.h
        src/double_node_object.h
        src/InputReader.cpp
        src/InputReader.h
        src/main.cpp
        src/node.cpp
        src/node.h
        src/node_object.h
        src/OptimizeBioGeo.cpp
        src/OptimizeBioGeo.h
        src/OptimizeBioGeoAllDispersal.cpp
        src/OptimizeBioGeoAllDispersal.h
        src/OptimizeBioGeoAllDispersal_nlopt.cpp
        src/OptimizeBioGeoAllDispersal_nlopt.h
        src/RateMatrixUtils.cpp
        src/RateMatrixUtils.h
        src/RateModel.cpp
        src/RateModel.h
        src/string_node_object.h
        src/superdouble.cpp
        src/superdouble.h
        src/tree.cpp
        src/tree.h
        src/tree_reader.cpp
        src/tree_reader.h
        src/tree_utils.cpp
        src/tree_utils.h
        src/Utils.cpp
        src/Utils.h
        src/vector_node_object.h)

target_link_libraries(lagrange_cpp gsl armadillo nlopt)
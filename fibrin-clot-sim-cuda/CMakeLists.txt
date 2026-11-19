cmake_minimum_required(VERSION 3.24)
# GPU architecture: 120 = RTX 50-series. Override with -DCMAKE_CUDA_ARCHITECTURES=89 (RTX 40), 86 (RTX 30), ...
set(CMAKE_CUDA_ARCHITECTURES 120 CACHE STRING "CUDA architectures")
project(fibrin_clot_sim LANGUAGES CXX CUDA)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)  # -O3; no fast-math
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
find_package(CUDAToolkit 12 REQUIRED)

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(PUGI ${SRC}/src_ext_lib/pugixml/src)
# System_Builder.cpp passes Thrust containers between translation units; compile it
# with nvcc too so Thrust's architecture-dependent inline namespace matches.
set_source_files_properties(${SRC}/System_Builder.cpp PROPERTIES LANGUAGE CUDA)

add_executable(bend-model
    ${SRC}/main.cu
    ${SRC}/System.cu
    ${SRC}/Storage.cu
    ${SRC}/System_Builder.cpp
    ${SRC}/src_bucket_scheme/Bucket_Net.cu
    ${SRC}/src_bucket_scheme/Bucket_Plt.cu
    ${SRC}/src_fibrin/Advance_Positions_Fibrin.cu
    ${SRC}/src_fibrin/Link_Nodes.cu
    ${SRC}/src_fibrin/WLC_Force.cu
    ${SRC}/src_fibrin/Torsion_Force.cu
    ${SRC}/src_fibrin/Params_Calc.cu
    ${SRC}/src_plt/Advance_Positions_Plt.cu
    ${SRC}/src_plt/Plt_Vol_Exc_Force.cu
    ${SRC}/src_plt/Plt_Field_Plt_Force.cu
    ${SRC}/src_plt/Plt_Field_Node_Force.cu
    ${SRC}/src_plt/Plt_Arm_Plt_Force.cu
    ${SRC}/src_plt/Plt_Arm_Node_Force.cu
    ${SRC}/src_plt/Plt_Arm_Node_Force_Time.cu
    ${PUGI}/pugixml.cpp
)
target_include_directories(bend-model PRIVATE
    ${SRC}
    ${SRC}/src_bucket_scheme
    ${SRC}/src_fibrin
    ${SRC}/src_plt
    ${SRC}/src_ext_lib/glm
    ${PUGI}
)
target_link_libraries(bend-model PRIVATE CUDA::cudart CUDA::curand)
target_compile_options(bend-model PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)

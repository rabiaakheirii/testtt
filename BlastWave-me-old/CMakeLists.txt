cmake_minimum_required(VERSION 3.16)
project(qcd_collision_sim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

add_executable(sim
    src/main.cpp
    src/core/config_parser.cpp
    src/utils/math.cpp
    src/utils/bessel.cpp
    src/utils/random.cpp
    src/freezeout/resonance_decay.cpp
    src/freezeout/blastwave.cpp
    src/analysis/spectra.cpp
    src/analysis/critical_point.cpp
    src/analysis/cumulants.cpp
    src/freezeout/cooper_frye.cpp
    src/hydro/hydro.cpp
    src/hydro/eos.cpp
    src/percolation/cluster.cpp
    src/percolation/percolation.cpp
    src/freezeout/longitudinal.cpp 
)

target_include_directories(sim PRIVATE src)
target_link_libraries(sim PRIVATE OpenMP::OpenMP_CXX GSL::gsl GSL::gslcblas)
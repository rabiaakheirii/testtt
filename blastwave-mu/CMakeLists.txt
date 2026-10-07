cmake_minimum_required(VERSION 3.16)
project(blastwave_mu CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

add_library(bwmu_core STATIC
    src/utils/bessel.cpp
    src/utils/integrate.cpp
    src/utils/grid.cpp
    src/utils/random.cpp
    src/utils/minimizer.cpp
    src/particles/table.cpp
    src/particles/chem_potential.cpp
    src/blastwave/spectrum.cpp
    src/blastwave/feeddown.cpp
    src/io/phenix_reader.cpp
    src/io/alice_reader.cpp
    src/io/mpd_reader.cpp
    src/io/csv_writer.cpp
    src/fit/chi2.cpp
    src/fit/chi2_surface.cpp
    src/fit/ptbin_scan.cpp
    src/fit/mcmc.cpp
    src/analysis/dndy.cpp
    src/fit/joint.cpp
    src/fit/tbeta_penalty.cpp
    src/fit/minuit_fit.cpp
)
target_include_directories(bwmu_core PUBLIC src)
target_link_libraries(bwmu_core PUBLIC
    OpenMP::OpenMP_CXX
    GSL::gsl
    GSL::gslcblas
)

add_executable(bwmu src/main.cpp)
target_link_libraries(bwmu PRIVATE bwmu_core)

enable_testing()

add_executable(test_phase1 tests/test_phase1.cpp)
target_link_libraries(test_phase1 PRIVATE bwmu_core)
add_test(NAME phase1 COMMAND test_phase1)

add_executable(test_phase2 tests/test_phase2.cpp)
target_link_libraries(test_phase2 PRIVATE bwmu_core)
add_test(NAME phase2 COMMAND test_phase2)

add_executable(test_phase3 tests/test_phase3.cpp)
target_link_libraries(test_phase3 PRIVATE bwmu_core)
add_test(NAME phase3 COMMAND test_phase3)

add_executable(test_phase4 tests/test_phase4.cpp)
target_link_libraries(test_phase4 PRIVATE bwmu_core)
add_test(NAME phase4 COMMAND test_phase4)

add_executable(test_phase5 tests/test_phase5.cpp)
target_link_libraries(test_phase5 PRIVATE bwmu_core)
add_test(NAME phase5 COMMAND test_phase5)

add_executable(test_phase6 tests/test_phase6.cpp)
target_link_libraries(test_phase6 PRIVATE bwmu_core)
add_test(NAME phase6 COMMAND test_phase6)

add_executable(test_phase8 tests/test_phase8.cpp)
target_link_libraries(test_phase8 PRIVATE bwmu_core)
add_test(NAME phase8 COMMAND test_phase8)
set_tests_properties(phase8 PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_phase10 tests/test_phase10.cpp)
target_link_libraries(test_phase10 PRIVATE bwmu_core)
add_test(NAME phase10 COMMAND test_phase10)

add_executable(test_phase11 tests/test_phase11.cpp)
target_link_libraries(test_phase11 PRIVATE bwmu_core)
add_test(NAME phase11 COMMAND test_phase11)

add_executable(test_phase12 tests/test_phase12.cpp)
target_link_libraries(test_phase12 PRIVATE bwmu_core)
add_test(NAME phase12 COMMAND test_phase12)

add_executable(test_phase13 tests/test_phase13.cpp)
target_link_libraries(test_phase13 PRIVATE bwmu_core)
add_test(NAME phase13 COMMAND test_phase13)


# ------------------------------------------------------------------
# ROOT (snap installation)
# ------------------------------------------------------------------
set(ROOT_INC "/snap/root-framework/current/usr/local/include")
set(ROOT_LIB "/snap/root-framework/current/usr/local/lib")

target_include_directories(bwmu_core PUBLIC ${ROOT_INC})
target_link_libraries(bwmu_core PUBLIC
    ${ROOT_LIB}/libMinuit2.so
    ${ROOT_LIB}/libCore.so
)

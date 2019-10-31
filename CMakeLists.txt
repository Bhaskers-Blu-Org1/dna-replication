cmake_minimum_required(VERSION 3.6)
project(dna-replication)

set(CMAKE_CXX_STANDARD 11)
if (NOT DEFINED EIGEN3_INCLUDE_DIR)
    set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
endif ()

include(GNUInstallDirs)

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system random filesystem program_options)

set(SOURCE_FILES        
    src/SimulationObserver.h
    src/Simulation.h src/Simulation.cpp
    src/Chromosome.h src/Chromosome.cpp
    src/Origin.h src/Origin.cpp
    src/Particle.h src/Particle.cpp
    src/ParticleDiffusionBehavior.h
    src/ParticleActivationBehavior.h
    src/ParticleBindingBehavior.h
    src/misc/Vector3.h
    src/misc/RandomTools.h
    src/misc/DataManager.h src/misc/DataManager.cpp
    src/simulation_observers/MultiSimulationObserver.h src/simulation_observers/MultiSimulationObserver.cpp
    src/particle_behaviors/IsotropicParticleDiffusionBehavior.h src/particle_behaviors/IsotropicParticleDiffusionBehavior.cpp
    src/particle_behaviors/ProbabilisticParticleActivationBehavior.h src/particle_behaviors/ProbabilisticParticleActivationBehavior.cpp
    src/particle_behaviors/ProbabilisticParticleBindingBehavior.h src/particle_behaviors/ProbabilisticParticleBindingBehavior.cpp)

add_executable(dnarepl src/SimulationClient.cpp ${SOURCE_FILES})
target_compile_definitions(dnarepl INTERFACE EIGEN_MPL2_ONLY)
target_link_libraries(dnarepl ${Boost_LIBRARIES})
target_include_directories(dnarepl PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
install(TARGETS dnarepl RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(dnarepl-batch src/SimulationBatchClient.cpp ${SOURCE_FILES})
target_compile_definitions(dnarepl-batch INTERFACE EIGEN_MPL2_ONLY)
target_compile_options(dnarepl-batch PUBLIC ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(dnarepl-batch ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(dnarepl-batch PUBLIC ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})
install(TARGETS dnarepl-batch RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
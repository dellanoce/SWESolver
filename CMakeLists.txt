cmake_minimum_required(VERSION 3.10)
project(SWE)

# -------------------------------------------- Compiler options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# ------------------------------------------------- Setup paths
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BIN_DIR ${ROOT_DIR}/bin)
set(SRC_DIR ${ROOT_DIR}/src)
set(INC_DIR ${ROOT_DIR}/include)

# -------------------------------------------------Source files
file(GLOB SOURCE ${SRC_DIR}/*.cpp
        ${INC_DIR}/*.h
        ${INC_DIR}/*.inl)

include_directories(${INC_DIR})

# ---------------------------------------------------- Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# --------------------------------------------------------- SWE
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
add_executable(${PROJECT_NAME} ${SOURCE})
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
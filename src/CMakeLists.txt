cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(yrt-pet-petsird VERSION 0.7.2)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(OpenMP REQUIRED)

#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif ()

if (WIN32)
    add_compile_options(/W3 /WX)
else ()
    add_compile_options(-Wall -Wextra -pedantic)
endif ()

set(CURRENT_PATH ./)
get_filename_component(CURRENT_PATH_REAL ${CURRENT_PATH} REALPATH)
set(PETSIRD_dir ${CURRENT_PATH_REAL}/../PETSIRD/cpp/generated)
get_filename_component(PETSIRD_dir_REAL ${PETSIRD_dir} REALPATH)
add_subdirectory(${PETSIRD_dir_REAL} PETSIRD_generated)

set(YRTPET_PETSIRD_SOURCES utils.cpp PETSIRDListMode.cpp DetectorCorrespondenceMap.cpp)

add_executable(petsird_yrtpet_reconstruct petsird_yrtpet_reconstruct.cpp ${YRTPET_PETSIRD_SOURCES})

target_link_libraries(petsird_yrtpet_reconstruct PUBLIC petsird_generated)
target_link_libraries(petsird_yrtpet_reconstruct PUBLIC yrtpet)
target_link_libraries(petsird_yrtpet_reconstruct PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(petsird_yrtpet_reconstruct PUBLIC ${PETSIRD_dir})
target_include_directories(petsird_yrtpet_reconstruct PUBLIC ${PETSIRD_dir}/..)
target_include_directories(petsird_yrtpet_reconstruct PUBLIC ${PETSIRD_dir}/../helpers/include)

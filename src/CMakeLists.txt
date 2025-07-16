cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(yrt-pet-petsird VERSION 0.7.2)

set(CMAKE_CXX_STANDARD 17)

#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif ()

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CURRENT_PATH ./)
get_filename_component(CURRENT_PATH_REAL ${CURRENT_PATH} REALPATH)
set(PETSIRD_dir ${CURRENT_PATH_REAL}/../PETSIRD/cpp/generated)
get_filename_component(PETSIRD_dir_REAL ${PETSIRD_dir} REALPATH)
add_subdirectory(${PETSIRD_dir_REAL} PETSIRD_generated)

add_executable(yrtpet_petsird_generate_scanner petsird_yrtpet_reconstruct.cpp)

target_link_libraries(yrtpet_petsird_generate_scanner PUBLIC petsird_generated)
target_include_directories(yrtpet_petsird_generate_scanner PUBLIC ${PETSIRD_dir})
target_include_directories(yrtpet_petsird_generate_scanner PUBLIC ${PETSIRD_dir}/..)
target_include_directories(yrtpet_petsird_generate_scanner PUBLIC ${PETSIRD_dir}/../helpers/include)

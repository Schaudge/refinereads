cmake_minimum_required(VERSION 3.13)
project(refinereads)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -Wno-unused-function -Wno-unused-variable -Wno-sign-compare)

### add extra (third-party) include directories
# include_directories(/opt/lib/include)
### add extra (third-party) lib directories
# link_directories(/opt/lib)

### add all source file to the SRC_LIST *variable*
aux_source_directory(src SRC_LIST)

add_executable(refinereads src/main.cpp ${SRC_LIST})
target_link_libraries(refinereads hts z pthread)

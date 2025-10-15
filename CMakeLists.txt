cmake_minimum_required(VERSION 3.22)
project(hausdorff)

set(CMAKE_CXX_STANDARD 14)

# ==== LIBRERÍA libcds ====
# Toma automáticamente todos los .cpp de libcds
file(GLOB_RECURSE LIBCDS_SRC
        ${CMAKE_SOURCE_DIR}/libcds/src/*.cpp
)

# Crea una librería estática llamada "cds"
add_library(cds STATIC ${LIBCDS_SRC})

# Rutas de include de libcds
target_include_directories(cds PUBLIC
        ${CMAKE_SOURCE_DIR}/libcds/include
        ${CMAKE_SOURCE_DIR}/libcds/build-aux
)

# ==== EJECUTABLE hausdorff ====
add_executable(hausdorff
        K2tree/K2tree.cpp
        Util/Factory.cpp
        Util/Point.cpp
        Util/TimeMesure.cpp
        utils.cpp
        ryu-kamata.cpp
        main.cpp
)

# Rutas de tus propios headers
target_include_directories(hausdorff PUBLIC
        .
        IO
        Logs
        K2tree
        Util
)

# Vincula tu ejecutable con la librería cds
target_link_libraries(hausdorff PRIVATE cds)
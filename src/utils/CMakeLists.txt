find_package(OpenMP)

add_library(PrivLR_UTILS
    data.cpp
    io.cpp
    paillier.cpp
)

target_include_directories(PrivLR_UTILS PRIVATE
    ${PROJECT_SOURCE_DIR}/include/utils/
)

target_link_libraries(PrivLR_UTILS PUBLIC
    ntl
    gmp
    OpenMP::OpenMP_CXX
    PrivLR
)
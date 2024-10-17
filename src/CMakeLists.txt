add_subdirectory(protocols)
add_subdirectory(utils)

add_library(PrivLR SHARED
    logistic-bfv.cpp
    logistic-paillier.cpp
)

target_include_directories(PrivLR PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(PrivLR PUBLIC
    PrivLR_PROTOCOLS
)
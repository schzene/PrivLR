add_library(PrivLR_PROTOCOLS
    linear-bfv.cpp
    linear-elgamal.cpp
    linear-paillier.cpp
    non-linear-bfv.cpp
    non-linear-paillier.cpp
)

target_include_directories(PrivLR_PROTOCOLS PRIVATE
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_SOURCE_DIR}/include/protocols
)

target_link_libraries(PrivLR_PROTOCOLS PUBLIC
    PrivLR_UTILS
)
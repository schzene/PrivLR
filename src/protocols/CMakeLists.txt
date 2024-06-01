add_library(PrivLR_PROTOCOLS
    linear.cpp
    non-linear.cpp
)

target_include_directories(PrivLR_PROTOCOLS PRIVATE
    ${PROJECT_SOURCE_DIR}/include/protocols
)

target_link_libraries(PrivLR_PROTOCOLS PUBLIC
    PrivLR_UTILS
)
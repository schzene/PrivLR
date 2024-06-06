include_directories(${PROJECT_SOURCE_DIR}/src)

macro(add_test _name)
    add_executable(${_name}_test "${_name}_test.cpp")
    target_link_libraries(${_name}_test PUBLIC PrivLR)
endmacro()

add_test(fixed-point)
add_test(linear)
add_test(logistic)
add_test(non-linear)
add_test(paillier)
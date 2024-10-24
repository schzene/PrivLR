include_directories(${PROJECT_SOURCE_DIR}/src)

macro(add_test _name)
    add_executable(${_name} "${_name}_test.cpp")
    target_link_libraries(${_name} PUBLIC privlr)
endmacro()

add_test(linear-bfv)
add_test(linear-elgamal)
add_test(linear-paillier)

add_test(logistic-bfv)
add_test(logistic-paillier)

add_test(non-linear-bfv)
add_test(non-linear-paillier)
include_directories(${PROJECT_SOURCE_DIR}/src)

macro(add_test _name)
    add_executable(${_name} "${_name}_test.cpp")
    target_link_libraries(${_name} PUBLIC PrivLR ntl)
endmacro()

add_test(logistic-bfv)
add_test(logistic-paillier)
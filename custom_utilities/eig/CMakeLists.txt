CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

add_library(eig STATIC eig3.cpp)
set_target_properties(eig PROPERTIES COMPILE_DEFINITIONS "KRATOS_CORE=IMPORT;STRUCTURAL_APPLICATION=EXPORT;")

set(file_list
test
)

foreach(str ${file_list})
    add_executable(${str} ${str}.cpp)
    target_link_libraries(${str} eig)
    set_target_properties(${str} PROPERTIES COMPILE_DEFINITIONS "KRATOS_CORE=IMPORT;STRUCTURAL_APPLICATION=NONE;")
endforeach()

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

include_directories(${CMAKE_SOURCE_DIR}/kratos)

add_library(eig STATIC eig3.cpp)
target_compile_definitions(eig PRIVATE KRATOS_CORE=IMPORT)
target_compile_definitions(eig PRIVATE STRUCTURAL_APPLICATION=EXPORT)

set(file_list
    test
)

foreach(str ${file_list})
    add_executable(${str} ${str}.cpp)
    target_link_libraries(${str} PRIVATE eig)
    target_compile_definitions(${str} PRIVATE KRATOS_CORE=IMPORT)
    target_compile_definitions(${str} PRIVATE STRUCTURAL_APPLICATION=NONE)
endforeach()

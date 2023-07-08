#message("Python version major: " ${PYTHON_VERSION_MAJOR})

message("")
message("**** configuring KratosStructuralApplication ****")

if(${PYTHON_VERSION_MAJOR} MATCHES "2")
    message("KratosStructuralApplication is compiled with Python 2")
    include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists2.txt)
elseif(${PYTHON_VERSION_MAJOR} MATCHES "3")
    message("KratosStructuralApplication is compiled with Python 3")
    include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists3.txt)
else()
    message("KratosStructuralApplication does not detect the Python version")
endif()

message("**** finish configuring KratosStructuralApplication ****")
message("")

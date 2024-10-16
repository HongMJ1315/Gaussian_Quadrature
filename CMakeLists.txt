cmake_minimum_required(VERSION 3.12)

set(MY_PROJECT "Hw1")
set(MY_EXECUTABLE "Hw1")

project("Hw1"
    VERSION 0.1.0
    LANGUAGES CXX C
    DESCRIPTION "Hw1"
)

add_executable(${MY_EXECUTABLE})
set_target_properties(${MY_EXECUTABLE}
    PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    C_STANDARD 11
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS OFF
)

# Find external libraries
find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(glm REQUIRED)
find_package(imgui REQUIRED)
find_package(OpenCV REQUIRED)

include_directories("include\\")
target_include_directories(${MY_EXECUTABLE} PRIVATE "include" ${STB_INCLUDE_DIRS})
file(GLOB MY_SOURCE CONFIGURE_DEPENDS
    "src/*.cpp"
    "src/*.cc"
    "src/*.c" # Add .c files to be compiled
)
target_sources(${MY_EXECUTABLE} PRIVATE ${MY_SOURCE})

target_link_libraries(${MY_EXECUTABLE} PRIVATE
    glfw
    GLEW::GLEW
    OpenGL::GL
    OpenGL::GLU
    GLUT::GLUT
    glad::glad
    glm::glm
    imgui::imgui
    opencv_ml
    opencv_dnn
    opencv_core
    opencv_flann
    opencv_imgcodecs
)

add_custom_command(TARGET ${MY_EXECUTABLE} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
    "${CMAKE_CURRENT_SOURCE_DIR}/shader"
    "$<TARGET_FILE_DIR:${MY_EXECUTABLE}>/shader"
    VERBATIM
)
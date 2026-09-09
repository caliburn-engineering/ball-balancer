# projects/ball-balancer/CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(ball-balancer LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default to an optimized build.  Measured (#10): a 200-point root-locus sweep
# at 10 closed-loop states costs 321 ms as this repo built before this line,
# versus 2.4 ms at -O2 -DNDEBUG — and needs_recompute is set from 26 sites
# including every slider-drag frame.  This makes NDEBUG live, which is only
# safe once the n == 0 guards in pole_zero/frequency_response exist, so the
# two land together.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)

# --- Shared reference implementations ---
# Golden-source C++ living in the Caliburn workspace, not in this repo.  Every
# project sits at the same depth under projects/, which is why ball-balancer's
# sources moved here without touching a single relative include.
set(CALIBURN_REFERENCE ${CMAKE_CURRENT_SOURCE_DIR}/../../reference)

# --- Eigen ---
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)

# --- GLFW (desktop only — Emscripten provides its own via -sUSE_GLFW=3) ---
if(NOT EMSCRIPTEN)
    FetchContent_Declare(glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG        3.4
        GIT_SHALLOW    TRUE
    )
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(glfw)
endif()

# --- GLAD (desktop only — Emscripten provides GL headers) ---
if(NOT EMSCRIPTEN)
    add_library(glad STATIC vendor/glad/src/glad.c)
    target_include_directories(glad PUBLIC vendor/glad/include)
endif()

# --- ImGui ---
FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        docking
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(imgui)
add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
if(EMSCRIPTEN)
    # USE_GLFW=3 provides GLFW headers via Emscripten ports (compile-time warning is harmless)
    target_compile_options(imgui PUBLIC -sUSE_GLFW=3)
else()
    target_link_libraries(imgui PUBLIC glfw glad)
    target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD)
endif()

# --- ImPlot ---
FetchContent_Declare(implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(implot)
add_library(implot STATIC
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
)
target_include_directories(implot PUBLIC ${implot_SOURCE_DIR})
target_link_libraries(implot PUBLIC imgui)

# --- Analysis library ---
add_library(analysis_lib STATIC
    ${CALIBURN_REFERENCE}/integrators/rk4.cpp
    src/analysis/model_library.cpp
    src/analysis/system_properties.cpp
    src/analysis/frequency_response.cpp
    src/analysis/pole_zero.cpp
    src/analysis/time_response.cpp
    src/analysis/system_connect.cpp
    src/analysis/controller_builders.cpp
    src/analysis/loop_diagnostics.cpp
    src/analysis/loop_locus.cpp
    src/analysis/lqr.cpp
)
target_include_directories(analysis_lib PUBLIC
    src
    ${CALIBURN_REFERENCE}/models       # linear_system.h, linearizer.h
    ${CALIBURN_REFERENCE}/integrators  # rk4.h
)
target_link_libraries(analysis_lib PUBLIC Eigen3::Eigen kinematics)

# --- Plate libraries ---
# Both targets, since issue #15: the kinematics, the rolling dynamics and the
# linearised plant are Eigen and the standard library only, so nothing here
# needed porting — only the renderer's shaders did.
add_library(kinematics STATIC
    src/table_kinematics.cpp
    src/trajectory_validator.cpp
)
target_include_directories(kinematics PUBLIC src)
target_link_libraries(kinematics PUBLIC Eigen3::Eigen)

add_library(ball_dynamics STATIC
    src/attract_mode.cpp
    src/ball_contact.cpp
    src/setpoint_path.cpp
    src/sim_step.cpp
    src/auto_balance.cpp
    src/ball_sim.cpp
    src/ball_plant_linear.cpp
    src/comparison_sim.cpp
    ${CALIBURN_REFERENCE}/simulation/rolling_dynamics.cpp
    ${CALIBURN_REFERENCE}/models/linearizer.cpp
)
target_include_directories(ball_dynamics PUBLIC
    src
    ${CALIBURN_REFERENCE}/simulation
)
# analysis_lib carries rk4 and LinearSystem for both halves of the app.
target_link_libraries(ball_dynamics PUBLIC analysis_lib kinematics)

# --- Tests (desktop only) ---
if(NOT EMSCRIPTEN)
    enable_testing()

    add_executable(test_model_library tests/test_model_library.cpp)
    target_link_libraries(test_model_library PRIVATE analysis_lib)
    target_include_directories(test_model_library PRIVATE tests)
    add_test(NAME test_model_library COMMAND test_model_library)

    add_executable(test_system_properties tests/test_system_properties.cpp)
    target_link_libraries(test_system_properties PRIVATE analysis_lib)
    target_include_directories(test_system_properties PRIVATE tests)
    add_test(NAME test_system_properties COMMAND test_system_properties)

    add_executable(test_frequency_response tests/test_frequency_response.cpp)
    target_link_libraries(test_frequency_response PRIVATE analysis_lib)
    target_include_directories(test_frequency_response PRIVATE tests)
    add_test(NAME test_frequency_response COMMAND test_frequency_response)

    add_executable(test_pole_zero tests/test_pole_zero.cpp)
    target_link_libraries(test_pole_zero PRIVATE analysis_lib)
    target_include_directories(test_pole_zero PRIVATE tests)
    add_test(NAME test_pole_zero COMMAND test_pole_zero)

    add_executable(test_time_response tests/test_time_response.cpp)
    target_link_libraries(test_time_response PRIVATE analysis_lib)
    target_include_directories(test_time_response PRIVATE tests)
    add_test(NAME test_time_response COMMAND test_time_response)

    add_executable(test_system_connect tests/test_system_connect.cpp)
    target_link_libraries(test_system_connect PRIVATE analysis_lib)
    target_include_directories(test_system_connect PRIVATE tests)
    add_test(NAME test_system_connect COMMAND test_system_connect)

    add_executable(test_controller_builders tests/test_controller_builders.cpp)
    target_link_libraries(test_controller_builders PRIVATE analysis_lib)
    target_include_directories(test_controller_builders PRIVATE tests)
    add_test(NAME test_controller_builders COMMAND test_controller_builders)

    add_executable(test_loop_diagnostics tests/test_loop_diagnostics.cpp)
    target_link_libraries(test_loop_diagnostics PRIVATE analysis_lib)
    target_include_directories(test_loop_diagnostics PRIVATE tests)
    add_test(NAME test_loop_diagnostics COMMAND test_loop_diagnostics)

    add_executable(test_lqr tests/test_lqr.cpp)
    target_link_libraries(test_lqr PRIVATE analysis_lib)
    target_include_directories(test_lqr PRIVATE tests)
    add_test(NAME test_lqr COMMAND test_lqr)

    add_executable(test_cascade_model tests/test_cascade_model.cpp)
    target_link_libraries(test_cascade_model PRIVATE analysis_lib)
    target_include_directories(test_cascade_model PRIVATE tests)
    add_test(NAME test_cascade_model COMMAND test_cascade_model)

    add_executable(test_ball_sim tests/test_ball_sim.cpp)
    target_link_libraries(test_ball_sim PRIVATE ball_dynamics)
    target_include_directories(test_ball_sim PRIVATE tests)
    add_test(NAME test_ball_sim COMMAND test_ball_sim)

    add_executable(test_auto_balance tests/test_auto_balance.cpp)
    target_link_libraries(test_auto_balance PRIVATE ball_dynamics)
    target_include_directories(test_auto_balance PRIVATE tests)
    add_test(NAME test_auto_balance COMMAND test_auto_balance)

    add_executable(test_trajectory tests/test_trajectory.cpp)
    target_link_libraries(test_trajectory PRIVATE ball_dynamics)
    target_include_directories(test_trajectory PRIVATE tests)
    add_test(NAME test_trajectory COMMAND test_trajectory)

    add_executable(test_setpoint_path tests/test_setpoint_path.cpp)
    target_link_libraries(test_setpoint_path PRIVATE ball_dynamics)
    target_include_directories(test_setpoint_path PRIVATE tests)
    add_test(NAME test_setpoint_path COMMAND test_setpoint_path)

    add_executable(test_ball_contact tests/test_ball_contact.cpp)
    target_link_libraries(test_ball_contact PRIVATE ball_dynamics)
    target_include_directories(test_ball_contact PRIVATE tests)
    add_test(NAME test_ball_contact COMMAND test_ball_contact)

    add_executable(test_assembly_mode tests/test_assembly_mode.cpp)
    target_link_libraries(test_assembly_mode PRIVATE ball_dynamics)
    target_include_directories(test_assembly_mode PRIVATE tests)
    add_test(NAME test_assembly_mode COMMAND test_assembly_mode)

    add_executable(test_attract_mode tests/test_attract_mode.cpp)
    target_link_libraries(test_attract_mode PRIVATE ball_dynamics)
    target_include_directories(test_attract_mode PRIVATE tests)
    add_test(NAME test_attract_mode COMMAND test_attract_mode)
endif()

# --- Visualizer executable ---
add_executable(visualizer
    src/visualizer.cpp
    src/panels/model_panel.cpp
    src/panels/properties_panel.cpp
    src/panels/pole_zero_panel.cpp
    src/panels/bode_panel.cpp
    src/panels/nyquist_panel.cpp
    src/panels/time_response_panel.cpp
    src/plate_view.cpp
    src/renderer.cpp
    src/plot_panel.cpp
    src/comparison_panel.cpp
)

if(EMSCRIPTEN)
    # GLFW, the GL headers and the ES 3.0 entry points all come from Emscripten
    # ports, so the web link needs no glad and no OpenGL::GL.
    target_link_libraries(visualizer PRIVATE
        analysis_lib
        ball_dynamics
        kinematics
        imgui
        implot
    )
    target_link_options(visualizer PRIVATE
        -sUSE_GLFW=3
        -sUSE_WEBGL2=1
        -sFULL_ES3=1
        -sALLOW_MEMORY_GROWTH=1
        -sMAX_WEBGL_VERSION=2
        --preload-file ${CMAKE_SOURCE_DIR}/vendor/fonts/NotoSans-Regular.ttf@vendor/fonts/NotoSans-Regular.ttf
        --shell-file ${CMAKE_SOURCE_DIR}/web/shell.html
    )
    set_target_properties(visualizer PROPERTIES SUFFIX ".html")
else()
    find_package(OpenGL REQUIRED)
    target_link_libraries(visualizer PRIVATE
        analysis_lib
        ball_dynamics
        kinematics
        imgui
        implot
        glad
        glfw
        OpenGL::GL
    )
endif()

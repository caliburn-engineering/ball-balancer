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
    src/rk4.cpp
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
target_include_directories(analysis_lib PUBLIC src)
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
    src/trajectory_controls.cpp
    src/sim_step.cpp
    src/loss_cause.cpp
    src/auto_balance.cpp
    src/ball_sim.cpp
    src/ball_plant_linear.cpp
    src/comparison_sim.cpp
    src/rolling_dynamics.cpp
    src/linearizer.cpp
)
target_include_directories(ball_dynamics PUBLIC src)
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

    add_executable(test_loss_cause tests/test_loss_cause.cpp)
    target_link_libraries(test_loss_cause PRIVATE ball_dynamics)
    target_include_directories(test_loss_cause PRIVATE tests)
    add_test(NAME test_loss_cause COMMAND test_loss_cause)

    add_executable(test_attract_mode tests/test_attract_mode.cpp)
    target_link_libraries(test_attract_mode PRIVATE ball_dynamics)
    target_include_directories(test_attract_mode PRIVATE tests)
    add_test(NAME test_attract_mode COMMAND test_attract_mode)

    # Wide pumping sweep — 21 gain nudges (1e-15 to 1e-12, both signs) over
    # the full 12x12 grid.  Off by default; run with
    #   cmake -DBALL_BALANCER_PUMPING_WIDE_SWEEP=ON
    # and then
    #   ctest -L slow
    # to re-measure the ceilings in test_the_loop_never_pumps_a_bouncing_ball.
    option(BALL_BALANCER_PUMPING_WIDE_SWEEP
           "Enable the 21-nudge wide sweep in test_attract_mode (slow, opt-in)"
           OFF)
    if(BALL_BALANCER_PUMPING_WIDE_SWEEP)
        target_compile_definitions(test_attract_mode
            PRIVATE BALL_BALANCER_PUMPING_WIDE_SWEEP)
        set_tests_properties(test_attract_mode PROPERTIES LABELS "slow")
    endif()

    # plate_design.cpp is compiled in so the guardrail test calls the real
    # designToOffer.  AppState names ImGui colour types, so the ImGui headers
    # are needed; the library is not linked, keeping GL out of the test.
    add_executable(test_guardrail tests/test_guardrail.cpp src/plate_design.cpp)
    target_link_libraries(test_guardrail PRIVATE ball_dynamics)
    target_include_directories(test_guardrail PRIVATE tests ${imgui_SOURCE_DIR})
    add_test(NAME test_guardrail COMMAND test_guardrail)

    add_executable(test_oracle_fixtures tests/test_oracle_fixtures.cpp)
    target_link_libraries(test_oracle_fixtures PRIVATE ball_dynamics)
    target_include_directories(test_oracle_fixtures PRIVATE tests tools)
    add_test(NAME test_oracle_fixtures COMMAND test_oracle_fixtures)
endif()

# --- Oracle exporter (desktop only, no GL needed) ---
# Runs the preset plants through the cascade model and emits A, B, Q, R as JSON
# so the Python oracle generator can solve CARE independently:
#   cmake --build build --target export_plants
#   ./build/export_plants | python3 tools/gen_oracle_fixtures.py - > tests/oracle_fixture.h
if(NOT EMSCRIPTEN)
    add_executable(export_plants tools/export_plants.cpp)
    target_link_libraries(export_plants PRIVATE ball_dynamics)
    target_include_directories(export_plants PRIVATE tools)
endif()

# --- Live-oracle option (desktop only) ---
# Regenerates oracle_fixture.h via the pipeline and cross-checks it against
# the checked-in version.  Requires python3 + scipy on the build machine;
# the test detects their absence at runtime and reports a ctest skip rather
# than a failure.
#
# Named with the top-level project prefix so the option survives promotion out
# of ball-balancer without a rename.
#
# Usage:
#   cmake -S . -B build -DCALIBURN_LIVE_ORACLE=ON \
#         -DCALIBURN_ORACLE_PYTHON=.oracle_venv/bin/python3
#   cmake --build build -j2
#   ctest --test-dir build -L live-oracle --output-on-failure
#
# CALIBURN_ORACLE_PYTHON is the interpreter the test runs the generator with;
# point it at the venv tools/requirements.txt was installed into.
if(NOT EMSCRIPTEN)
    option(CALIBURN_LIVE_ORACLE
           "Enable live-oracle tests: regenerate + cross-check oracle_fixture.h (requires python3 + scipy)"
           OFF)
    set(CALIBURN_ORACLE_PYTHON "python3" CACHE STRING
        "Python interpreter with scipy, used by the live-oracle tests")
    if(CALIBURN_LIVE_ORACLE)
        add_executable(test_live_oracle tests/test_live_oracle.cpp)
        target_compile_definitions(test_live_oracle PRIVATE
            EXPORT_PLANTS_EXE="$<TARGET_FILE:export_plants>"
            GEN_ORACLE_SCRIPT="${CMAKE_SOURCE_DIR}/tools/gen_oracle_fixtures.py"
            ORACLE_FIXTURE_FILE="${CMAKE_SOURCE_DIR}/tests/oracle_fixture.h"
            ORACLE_PYTHON="${CALIBURN_ORACLE_PYTHON}"
        )
        add_dependencies(test_live_oracle export_plants)
        add_test(NAME test_live_oracle COMMAND test_live_oracle)
        set_tests_properties(test_live_oracle PROPERTIES
            LABELS       "live-oracle"
            SKIP_RETURN_CODE 77)

        # Second independent oracle: pydrake's ContinuousAlgebraicRiccatiEquation
        # (Schur-based / LAPACK).  Skipped at runtime when pydrake is absent;
        # no new build option — it rides under CALIBURN_LIVE_ORACLE.
        add_executable(test_pydrake_oracle tests/test_pydrake_oracle.cpp)
        target_compile_definitions(test_pydrake_oracle PRIVATE
            EXPORT_PLANTS_EXE="$<TARGET_FILE:export_plants>"
            PYDRAKE_ORACLE_SCRIPT="${CMAKE_SOURCE_DIR}/tools/gen_oracle_pydrake.py"
            ORACLE_FIXTURE_FILE="${CMAKE_SOURCE_DIR}/tests/oracle_fixture.h"
            ORACLE_PYTHON="${CALIBURN_ORACLE_PYTHON}"
        )
        add_dependencies(test_pydrake_oracle export_plants)
        add_test(NAME test_pydrake_oracle COMMAND test_pydrake_oracle)
        set_tests_properties(test_pydrake_oracle PROPERTIES
            LABELS       "live-oracle"
            SKIP_RETURN_CODE 77)
    endif()
endif()

# --- Visualizer executable ---
add_executable(visualizer
    src/visualizer.cpp
    src/plate_design.cpp
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

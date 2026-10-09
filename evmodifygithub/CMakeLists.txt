cmake_minimum_required(VERSION 3.10)
project(ev_modify_solver VERSION 1.0.0 LANGUAGES C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wdouble-promotion -Werror)
endif()

# Solver library (shared, for external tools)
add_library(ev_modify_solver_shared SHARED spme_v7_apex_solver.c)

# Demo simulator + self-test
add_executable(ev_modify_solver spme_v7_apex_solver.c test_solver.c)

# stdin/stdout OCV proxy used by run_pybamm_sim.py and test_ocv.py
add_executable(ocv_proxy spme_v7_apex_solver.c ocv_proxy.c)

if(NOT MSVC)
    target_link_libraries(ev_modify_solver_shared PRIVATE m)
    target_link_libraries(ev_modify_solver PRIVATE m)
    target_link_libraries(ocv_proxy PRIVATE m)
endif()

enable_testing()
add_test(NAME CControllerSelfTest COMMAND ev_modify_solver --self-test)

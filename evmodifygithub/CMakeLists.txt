cmake_minimum_required(VERSION 3.10)
project(ev_modify_solver VERSION 1.0.0 LANGUAGES C)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
add_executable(ev_modify_solver spme_v7_apex_solver.c)
add_library(ev_modify_solver_shared SHARED spme_v7_apex_solver.c)
target_compile_definitions(ev_modify_solver_shared PRIVATE EV_TEST_MODE) # To expose STATIC functions if needed
if(NOT MSVC)
    target_link_libraries(ev_modify_solver PRIVATE m)
    target_link_libraries(ev_modify_solver_shared PRIVATE m)
endif()

enable_testing()
add_test(NAME CControllerSelfTest COMMAND ev_modify_solver --self-test)

include(FetchContent)
FetchContent_Declare(
  unity
  GIT_REPOSITORY https://github.com/ThrowTheSwitch/Unity.git
  GIT_TAG        v2.5.2
)
FetchContent_MakeAvailable(unity)

add_executable(test_solver_c test_solver_c.c spme_v7_apex_solver.c)
target_compile_definitions(test_solver_c PRIVATE EV_TEST_MODE)
target_link_libraries(test_solver_c PRIVATE unity)
if(NOT MSVC)
    target_link_libraries(test_solver_c PRIVATE m)
endif()
add_test(NAME CControllerUnityTest COMMAND test_solver_c)
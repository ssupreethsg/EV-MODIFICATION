# Minimal consumer of ev_modify, used by CI to prove the install package (find_package) and the
# add_subdirectory / FetchContent route both work.
#   find_package: cmake -S tests/consumer -B b -DCMAKE_PREFIX_PATH=<install prefix>
#   in-tree     : cmake -S tests/consumer -B b -DEV_CONSUMER_SOURCE_DIR=<repository root>
cmake_minimum_required(VERSION 3.21)
project(ev_modify_consumer LANGUAGES C)

set(EV_CONSUMER_SOURCE_DIR "" CACHE PATH "Repository root to add_subdirectory() instead of find_package()")
if(EV_CONSUMER_SOURCE_DIR)
    add_subdirectory(${EV_CONSUMER_SOURCE_DIR} ev_modify_build)
else()
    find_package(ev_modify 2.0 REQUIRED)
endif()

add_executable(consumer main.c)
target_link_libraries(consumer PRIVATE ev_modify::ev_modify)

enable_testing()
add_test(NAME consumer COMMAND consumer)

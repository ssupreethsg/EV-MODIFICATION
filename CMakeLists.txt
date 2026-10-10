cmake_minimum_required(VERSION 3.21)
project(ev_modify VERSION 2.0.0 LANGUAGES C)

include(GNUInstallDirs)

option(EV_BUILD_TESTS "Build the self-test, the CLI bridge and register the tests" ${PROJECT_IS_TOP_LEVEL})
option(EV_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (CI turns this on)" OFF)
option(EV_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer (GCC/Clang)" OFF)
option(EV_INSTALL "Generate the install and package rules" ${PROJECT_IS_TOP_LEVEL})
if(EV_INSTALL AND EV_SANITIZE)
    message(WARNING "EV_SANITIZE builds are for testing only: install rules disabled")
    set(EV_INSTALL OFF)
endif()
option(EV_PYTHON_TESTS "Register the Python checks with CTest (need numpy and scipy)" ON)

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)
    if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
endif()

if(MSVC)
    set(EV_WARNING_FLAGS /W4)
    if(EV_WARNINGS_AS_ERRORS)
        list(APPEND EV_WARNING_FLAGS /WX)
    endif()
else()
    set(EV_WARNING_FLAGS -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wdouble-promotion)
    if(EV_WARNINGS_AS_ERRORS)
        list(APPEND EV_WARNING_FLAGS -Werror)
    endif()
endif()

set(EV_SANITIZER_FLAGS)
if(EV_SANITIZE AND MSVC)
    message(WARNING "EV_SANITIZE is ignored with MSVC (only GCC and Clang are supported)")
endif()
if(EV_SANITIZE AND NOT MSVC)
    set(EV_SANITIZER_FLAGS -fsanitize=address,undefined,float-cast-overflow -fno-sanitize-recover=all -fno-omit-frame-pointer)
endif()

# ---- Controller library -------------------------------------------------------------------
add_library(ev_modify STATIC src/ev_modify.c)
add_library(ev_modify::ev_modify ALIAS ev_modify)
target_compile_features(ev_modify PUBLIC c_std_11)
target_include_directories(ev_modify PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_options(ev_modify PRIVATE ${EV_WARNING_FLAGS} ${EV_SANITIZER_FLAGS})
target_link_options(ev_modify INTERFACE "$<BUILD_INTERFACE:${EV_SANITIZER_FLAGS}>")
if(NOT MSVC)
    target_link_libraries(ev_modify PUBLIC m)
endif()

# ---- Self-test, CLI bridge and tests ------------------------------------------------------
if(EV_BUILD_TESTS)
    add_executable(ev_modify_selftest tests/test_solver.c)
    set_target_properties(ev_modify_selftest PROPERTIES OUTPUT_NAME test_solver)
    target_compile_options(ev_modify_selftest PRIVATE ${EV_WARNING_FLAGS} ${EV_SANITIZER_FLAGS})
    target_link_libraries(ev_modify_selftest PRIVATE ev_modify)

    add_executable(ev_modify_cli tools/ev_cli.c)
    set_target_properties(ev_modify_cli PROPERTIES OUTPUT_NAME ev_cli)
    target_compile_options(ev_modify_cli PRIVATE ${EV_WARNING_FLAGS} ${EV_SANITIZER_FLAGS})
    target_link_libraries(ev_modify_cli PRIVATE ev_modify)

    enable_testing()
    add_test(NAME ev_modify.c_self_test COMMAND ev_modify_selftest --self-test)
    add_test(NAME ev_modify.c_demo COMMAND ev_modify_selftest)

    # CLI contract: exit 2 = usage / argument error, exit 3 = malformed or out-of-range input line.
    set(EV_CLI_CELL "1.0|0.05|0.02|30")
    set(EV_CLI_PACK "${EV_CLI_CELL}|0.5|25")
    function(ev_add_cli_test name expect cli_args)
        add_test(NAME ev_modify.cli_${name}
                 COMMAND ${CMAKE_COMMAND} -DCLI=$<TARGET_FILE:ev_modify_cli> "-DCLI_ARGS=${cli_args}"
                         -DEXPECT=${expect} ${ARGN} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli/check_cli.cmake)
    endfunction()
    ev_add_cli_test(ok_ekf 0 "ekf|${EV_CLI_PACK}" "-DINPUT=1 0 3.7 25\\n")
    ev_add_cli_test(no_args 2 "")
    ev_add_cli_test(unknown_mode 2 "bogus")
    ev_add_cli_test(missing_args 2 "ekf|1.0")
    ev_add_cli_test(arg_not_a_number 2 "ekf|abc|0.05|0.02|30|0.5|25")
    ev_add_cli_test(arg_out_of_range 2 "ekf|1.0|0.05|0.02|30|1.5|25")
    ev_add_cli_test(arg_huge 2 "limit|1.0|0.05|0.02|30|1e300")
    ev_add_cli_test(bad_cell_config 2 "ekf|-1.0|0.05|0.02|30|0.5|25")
    ev_add_cli_test(ocv_not_a_number 3 "ocv" "-DINPUT=abc\\n")
    ev_add_cli_test(ocv_out_of_range 3 "ocv" "-DINPUT=1.5\\n")
    ev_add_cli_test(ekf_too_few_fields 3 "ekf|${EV_CLI_PACK}" "-DINPUT=1 0 3.7\\n")
    ev_add_cli_test(ekf_trailing_garbage 3 "ekf|${EV_CLI_PACK}" "-DINPUT=1 0 3.7 25 x\\n")
    ev_add_cli_test(ekf_dt_too_large 3 "ekf|${EV_CLI_PACK}" "-DINPUT=30 1 3.7 25\\n")
    ev_add_cli_test(ekf_dt_zero 3 "ekf|${EV_CLI_PACK}" "-DINPUT=0 1 3.7 25\\n")
    ev_add_cli_test(ekf_dt_too_small 3 "ekf|${EV_CLI_PACK}" "-DINPUT=0.0001 1 3.7 25\\n")
    ev_add_cli_test(ekf_voltage_huge 3 "ekf|${EV_CLI_PACK}" "-DINPUT=1 0 1e300 25\\n")
    ev_add_cli_test(charge_dt_too_large 3 "charge|${EV_CLI_PACK}" "-DINPUT=30 1 3.7 25 10 30\\n")
    ev_add_cli_test(limit_soc_out_of_range 3 "limit|${EV_CLI_CELL}|25" "-DINPUT=2 0.01 25 0 1 10 30\\n")
    ev_add_cli_test(nul_byte 3 "ocv" -DINPUT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/tests/cli/nul_line.bin)

    if(EV_PYTHON_TESTS)
        find_package(Python3 COMPONENTS Interpreter)
        set(EV_PYTHON_OK FALSE)
        if(Python3_FOUND)
            execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy"
                            RESULT_VARIABLE EV_PYTHON_IMPORT OUTPUT_QUIET ERROR_QUIET)
            if(EV_PYTHON_IMPORT EQUAL 0)
                set(EV_PYTHON_OK TRUE)
            endif()
        endif()
        if(EV_PYTHON_OK)
            add_test(NAME ev_modify.python_ocv
                     COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ocv.py --cli $<TARGET_FILE:ev_modify_cli>)
            add_test(NAME ev_modify.python_ocv_tables_current
                     COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_ocv_lut.py --check)
            add_test(NAME ev_modify.python_calibration
                     COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/tools/calibrate_cell.py --self-test)
            # keep __pycache__ out of the source tree
            set_tests_properties(ev_modify.python_ocv ev_modify.python_ocv_tables_current ev_modify.python_calibration
                                 PROPERTIES ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
        else()
            message(WARNING "Python 3 with numpy and scipy not found: Python checks not registered "
                            "(pip install -r requirements.txt, or -DEV_PYTHON_TESTS=OFF)")
        endif()
    endif()
endif()

# ---- Install / package --------------------------------------------------------------------
if(EV_INSTALL)
    include(CMakePackageConfigHelpers)
    install(TARGETS ev_modify EXPORT ev_modifyTargets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(FILES src/ev_modify.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(EXPORT ev_modifyTargets NAMESPACE ev_modify:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ev_modify)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ev_modifyConfig.cmake.in
         "@PACKAGE_INIT@\ninclude(\"\${CMAKE_CURRENT_LIST_DIR}/ev_modifyTargets.cmake\")\ncheck_required_components(ev_modify)\n")
    configure_package_config_file(${CMAKE_CURRENT_BINARY_DIR}/ev_modifyConfig.cmake.in
                                  ${CMAKE_CURRENT_BINARY_DIR}/ev_modifyConfig.cmake
                                  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ev_modify)
    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ev_modifyConfigVersion.cmake
                                     COMPATIBILITY SameMajorVersion)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ev_modifyConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/ev_modifyConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ev_modify)
    install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()

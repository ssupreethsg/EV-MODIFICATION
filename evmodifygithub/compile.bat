@echo off
set CFLAGS=-std=c11 -O2 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wdouble-promotion -Werror

echo Compiling EV MODIFY Physics Library...
gcc %CFLAGS% -shared -o spme_v7_apex_solver.dll spme_v7_apex_solver.c
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compile DLL
    exit /b %ERRORLEVEL%
)

echo Compiling C Simulator / Self-Test...
gcc %CFLAGS% -o ev_modify_solver.exe spme_v7_apex_solver.c test_solver.c
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compile test binary
    exit /b %ERRORLEVEL%
)

echo Compiling OCV proxy (stdin/stdout bridge used by the Python scripts)...
gcc %CFLAGS% -o ocv_proxy.exe spme_v7_apex_solver.c ocv_proxy.c
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compile proxy
    exit /b %ERRORLEVEL%
)

echo Running self-test...
ev_modify_solver.exe --self-test
if %ERRORLEVEL% NEQ 0 (
    exit /b %ERRORLEVEL%
)

echo Compilation Successful.

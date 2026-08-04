@echo off
echo Compiling EV MODIFY Physics Library...
gcc -shared -o spme_v7_apex_solver.dll -fPIC spme_v7_apex_solver.c
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compile DLL
    exit /b %ERRORLEVEL%
)

echo Compiling C Self-Test...
gcc -o test_bin.exe spme_v7_apex_solver.c test_solver.c
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compile test binary
    exit /b %ERRORLEVEL%
)

echo Compiling Python proxy (for 64-bit Python to 32-bit GCC interop)...
gcc -o ocv_proxy.exe spme_v7_apex_solver.c ocv_proxy.c
if %ERRORLEVEL% NEQ 0 (
    echo Failed to compile proxy
    exit /b %ERRORLEVEL%
)

echo Compilation Successful.

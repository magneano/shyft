@echo off
@rem simple build script that
@rem a) build shyft core debug and deploy it to SHYFT_DEPENDENCIES
@rem b) build shyft full release
@rem preconditions: shyft dependencies are built, and env. variables according to this is set(typical dev.setup)

pushd "%~dp0"\..
if defined VS2017INSTALLDIR (
    set msb="%VS2017INSTALLDIR%\MSbuild\15.0\bin\amd64\msbuild.exe"
) else (
    if defined MSBUILD_2017_PATH (
    set msb=%MSBUILD_2017_PATH%
    ) else (
    echo "ERROR NOT ABLE TO FIGURE OUT PATH TO MSBUILD, BAILING OUT"
    popd
    exit /b -1
    )
)
if not defined SHYFT_DEPENDENCIES (
    set SHYFT_DEPENDENCIES="%~dp0"\..\..\shyft_dependencies"
)
%msb% /m /p:Configuration=Debug /p:Platform=x64 /p:PlatformToolset=v141 /p:WindowsTargetPlatformVersion=10.0.16299.0 -target:core
if not %errorlevel%==0 (
    goto finale
) else (
    robocopy /NS /NFL /NP /NJH core %SHYFT_DEPENDENCIES%\include\shyft\core *.h
    copy /Y bin\Debug\shyft_core.lib %SHYFT_DEPENDENCIES%\lib\shyft_core.Debug.lib
    copy /Y bin\Debug\shyft_core.pdb %SHYFT_DEPENDENCIES%\lib\shyft_core.Debug.pdb
)

%msb% /m /p:Configuration=Release /p:Platform=x64 /p:PlatformToolset=v141 /p:WindowsTargetPlatformVersion=10.0.16299.0
if not %errorlevel%==0 (
    goto finale
) else (
    copy /Y bin\Release\shyft_core.lib %SHYFT_DEPENDENCIES%\lib\shyft_core.Release.lib
)
goto finale

:finale
popd
exit /b %errorlevel%
@echo off
@rem simple conda build/upload script that
@rem a) set version
@rem b) conda build (and upload if you have set upload to true)
@rem preconditions: shyft is built and tested

pushd "%~dp0"\..

git rev-list --count HEAD >VERSION.tmp
set /p minor=<VERSION.tmp
echo 4.6.%minor% >VERSION
set SHYFT_VERSION=4.6.%minor%
echo Starting to build %SHYFT_VERSION% numpy %SHYFT_BOOST_NUMPY_VERSION%
conda build --numpy %SHYFT_BOOST_NUMPY_VERSION% --no-test --no-copy-test-source-files conda_recipe
popd
exit /b %errorlevel%
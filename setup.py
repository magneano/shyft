import os
import re
import sys
import platform
import shutil
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

if sys.version_info < (3, 6):
    print("Python 3.6 or higher required, please upgrade.")
    sys.exit(1)

if os.path.exists('VERSION'):
    VERSION = open('VERSION').read().strip()
    # Create the version.py file
    open('shyft/version.py', 'w').write('__version__ = "%s"\n'%VERSION)
else:
    from shyft.version import __version__ as VERSION

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='',
                 files=None):
        Extension.__init__(self, name, sources=[])
        self.files = files if files is not None else {}
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                         out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def copy_files(self, ext):
        build_dir = os.path.realpath(self.build_temp)
        root_dir = os.path.dirname(os.path.realpath(__file__))
        lib_dir = os.path.realpath(os.path.join(self.build_lib, ext.name))

        for src_file, dest_file in ext.files.items():
            src = os.path.join(build_dir, src_file)
            dest = os.path.join(lib_dir, dest_file)

            if not os.path.exists(os.path.dirname(dest)):
                os.makedirs(os.path.dirname(dest))

            print(f"copying file {src} -> {dest}")
            shutil.copyfile(src, dest)

    def build_extension(self, ext):
        cmake_args = ['-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)

        self.copy_files(ext)
        print()  # Add an empty line for cleaner output

ext_s = 'pyd' if 'Windows' in platform.platform() else 'so'
shyft_files = files={f'api/boostpython/_api.{ext_s}': f'api/_api.{ext_s}',
                     f'api/boostpython/_hbv_stack.{ext_s}': f'api/hbv_stack/_hbv_stack.{ext_s}',
                     f'api/boostpython/_pt_gs_k.{ext_s}': f'api/pt_gs_k/_pt_gs_k.{ext_s}',
                     f'api/boostpython/_pt_hps_k.{ext_s}': f'api/pt_hps_k/_pt_hps_k.{ext_s}',
                     f'api/boostpython/_pt_hs_k.{ext_s}': f'api/pt_hs_k/_pt_hs_k.{ext_s}',
                     f'api/boostpython/_pt_ss_k.{ext_s}': f'api/pt_ss_k/_pt_ss_k.{ext_s}'}

setup(name='shyft',
      version=VERSION,
      author='Statkraft',
      description='An OpenSource hydrological toolbox',
      long_description='',
      url='https://github.com/statkraft/shyft',
      license='LGPL v3',
      packages=find_packages(),
      ext_modules=[CMakeExtension('shyft', files=shyft_files)],
      cmdclass=dict(build_ext=CMakeBuild),
      install_requires=['numpy'],
      tests_require=['nose'],
      extras_require={'repositories': ['netcdf4', 'shapely', 'pyyaml', 'six', 'pyproj'],
                      'notebooks': ['jupyter']},
      entry_points={},
      zip_safe=False)

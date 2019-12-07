import sys, os
from cx_Freeze import setup, Executable
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))

sys.path.append('pyamrfinder')

includes = ["pyamrfinder"]
includefiles = ["pyamrfinder/logo.gif","pyamrfinder/data",
                os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
                os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll')
                ]

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["os","numpy","matplotlib",
                                  "pandas","pandastable",
                                  "pyamrfinder"],
                     "excludes": ['scipy','seaborn','statsmodels'],
                     "namespace_packages": ['mpl_toolkits'],
                     "include_msvcr": True,
                     "includes": includes,
                     "include_files": includefiles}

base = None
if sys.platform == "win32":
    base = "Win32GUI"

executables = [Executable("main.py", base=base,
                          #copyDependentFiles = True,
                          targetName='pyamrfinder.exe',
                          shortcutName="pyamrfinder",
                          shortcutDir="pyamrfinder",
                          icon="img/logo.ico")]

setup(  name = "pyamrfinder",
	version = "0.1.0",
	description = "",
    options = {"build_exe": build_exe_options},
    executables = executables)

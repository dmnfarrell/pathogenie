import sys, os
from cx_Freeze import setup, Executable
PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))

sys.path.append('pygenefinder')

includes = ["pygenefinder"]
includefiles = ["pygenefinder/logo.gif","pygenefinder/data",
                os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
                os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll'),                
                ]

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["os","numpy","matplotlib","Bio",
                                  "pandas","pandastable",
                                  "pygenefinder"],
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
                          targetName='pygenefinder.exe',
                          shortcutName="pygenefinder",
                          shortcutDir="pygenefinder",
                          icon="img/logo.ico")]

setup(  name = "pygenefinder",
	version = "0.1.0",
	description = "amr gene finder",
    options = {"build_exe": build_exe_options},
    executables = executables)

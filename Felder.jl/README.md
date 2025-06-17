# <div align="center"><img src="/images/logo_96dpi_transparent.png?raw=true" width="350" title="Felder logo"></div>

Multiphysics FEM in 1D, 2D and 3D.

Experimental code. Use at your own risk.

Created and tested in Julia 1.9.4.

# 1. Julia Installation

1. Download the latest Julia installer from https://julialang.org/downloads/.
2. Run the installer and following the instructions on the screen.

Make sure that the Julia executable is added to the PATH environment variable. You can then
start the Julia interactive command-line REPL (read-eval-print loop) by executing `julia`
from any directory in a command-line (e.g. PowerShell, cmd, Bash, ...).

For more information, visit

https://docs.julialang.org/en/v1/stdlib/REPL/

# 2. Setting up the Felder.jl Julia environment (first execution only)

Before you can run the Felder.jl package for the first time, all necessary dependencies (specified in Project.toml)
need to be installed into a local Julia environment. This can be achieved in two ways (2.A. **or** 2.B.):

## 2.A. Directly from the command-line

1. Open the `.\Felder\` directory in a command-line prompt.
2. Enter `julia --project --eval 'using Pkg; Pkg.instantiate()'`.

## 2.B. Interactively from the Julia REPL

1. Open the `.\Felder\` directory in a command-line prompt.
3. Enter `julia --project` to start the Julia REPL in a local environment.
4. Enter `]instantiate` to invoke the package manager and install all dependencies.
5. Close the Julia session (`CTRL + D`).

For more information on environments, visit

https://pkgdocs.julialang.org/v1/environments/

# 3. Running a script

Felder.jl simulation scripts are located in the `.\Felder\scripts\` directory and are intended to be edited by the user.
There are two ways to run a script (3.A. **or** 3.B.):

## 3.A. Directly from the command-line

1. Open the `.\Felder\` directory in a command-line prompt.
2. Enter `julia --project --threads=auto .\scripts\myscript.jl` to run the script with multi-threading enabled.

## 3.B. Interactively from the Julia REPL

This is the recommended way when you are editing and exploring a script and plan to execute it multiple times. It keeps the Julia instance running so that the Felder source code only needs to be precompiled once.

1. Open the `.\Felder\` directory in a command-line prompt.
2. Enter `julia --project --threads=auto` to start the Julia REPL in the local environment with multi-threading enabled.
3. Enter `include("scripts/myscript.jl")` to run the script.

## 3.1. Multithreading

The above commands to run a script already contain the command-line argument `--threads=auto` that enables multi-threading. "auto" sets the number of threads to the number of local CPU threads.

Alternatively you can create an environment variable "JULIA_NUM_THREADS" which specifies the number of threads that will be used when the `--threads=` argument is not present.

For more details, visit

https://docs.julialang.org/en/v1/manual/multi-threading/

## 3.2. More command-line options

For more command-line options, visit

https://docs.julialang.org/en/v1/manual/command-line-options/

## 3.3 Julia REPL help mode

You can browse the Julia and Felder.jl documentation from the Julia REPL by entering the help mode by pressing `?` and then entering the name of the function or type you want to learn more about. For example, in a running REPL session:

1. Enter `using Felder` to load the Felder.jl package if necessary.
2. Enter `?LinearInterpolation` to print the `LinearInterpolation` documentation.

For more informatin, visit

https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode

# 4. Visualization of the solution in ParaView

Executing the `solve!()` function will write binary .vti files (Visualization Toolkit for Unstructured Grids) with a .pvd time series file into an output folder.

Each .vtu file contains the solution data (e.g. temperature, stress, ...) including the mesh data (e.g. node coordinates, element connectivity, ...) at a specific time step.

The .pvd file contains a list of all .vtu files and their corresponding time steps so that they can be visualized as a time series in ParaView.

To download ParaView, visit:

https://www.paraview.org/download/

## 4.A. Creating a custom visualization

1. Open ParaView and click "File" -> "Open...".
2. Select a .vti file to open a single time step or a .pvd file to open a time series.
3. In the "Pipeline Browser" on the left hand side, click on the eye symbol next to the opened file to render the data in the default view.
4. In the second tool row from above, open the drop down menu that shows "Solid Color" to select a solution variable to visualize.
5. In the first tool row from above you can find time controls to play the time series.

For a comprehensive introduction to ParaView and its visualization capabilities, visit

https://www.paraview.org/tutorials/

## 4.1. Saving a state file (XML or Python script)

ParaView visualzation settings (state) can be saved in a .pvsm file or a .py Python script. However, the .pvsm and the default .py format contain absolute file paths to the solution files, which will not work on other computers unless the Python scritps are edited accordingly or the .pvsm file is opened as described in section 4.C.

1. Open ParaView and click "File" -> "Save state...".
2. Choose a file name and type (.py or .pvsm) and click "Save".

# 5. Development Tips

If you are developing the source code, a "Revised-based" workflow is recommended to avoid the need for precompilation after every code change. This is supported by default in the Julia VSCode extension.

For more information, visit

https://docs.julialang.org/en/v1/manual/workflow-tips/#Revise-based-workflows

https://www.julia-vscode.org/


## Set Up Environment (to run in WSL)

### Install conda
*Note: These instructions are for installing conda in Linux.*

Run the following code in the Ubuntu terminal:

```bash
# Download the installer
# Replace <INSTALLER_VERSION> with the version you would like (I used 2024.02-1)
curl -O https://repo.anaconda.com/archive/Anaconda3-<INSTALLER_VERSION>-Linux-x86_64.sh

# Install
bash Anaconda3-<INSTALLER_VERSION>-Linux-x86_64.sh
```

### Create Environment
After installing conda, create the environment by running the following code in the terminal:

```bash
# Create the environment
# Replace <ENV_NAME> with a name of your choice (I used "heatequationenv")
# Replace <PYTHON_VERSION> with a version of your choice (I used 3.12.2)
conda create -n <ENV_NAME> python=<PYTHON_VERSION>

# Activate the environment
conda activate <ENV_NAME>
```

### Install FEniCS, matplotlib, and Jupyter Notebook
*Note: This installation method is only available for conda.*

1. After activating the environment, run in the terminal:

```bash
conda install -c conda-forge fenics
conda install -c conda-forge matplotlib
conda install jupyter
```

2. Close the terminal.

### Open Jupyter Notebook
To create a directory for the project and open Jupyter Notebook for writing code, open a terminal and run the following code:

```bash
# Make directory
# Replace <DIRECTORY_NAME> with a name of your choice
conda activate <ENV_NAME>
mkdir <DIRECTORY_NAME>

# Open Jupyter Notebook
jupyter notebook
```

Running `jupyter notebook` creates one HTML file path and two URLs. Copy and paste one of the URLs or HTML file path into a browser to access Jupyter Notebook (with WSL).


## Set Up Solutions for Solving the Heat Equation with the Finite Element Method
The entire code for solving the heat equation using the finite element method is:

```python
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Set parameters
k = 0.1
alpha = 3
nx, ny = 100,100

# Create mesh and define function space
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'CG', 2)

# Define boundary conditions
# Boundary conditions will be u(x,y) = 1+x^2+alpha*y^2
u_D = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1]', 
                degree=2, alpha=alpha)

def boundary(x, on_boundary):
    return on_boundary

# Set boundary conditions
bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
q = Constant(-k*(2+2*alpha))
n = FacetNormal(mesh)
F = k*dot(grad(u), grad(v))*dx - k*v*dot(grad(u),n)*ds - q*v*dx
a, L = lhs(F), rhs(F)

# Solve
u = Function(V)
solve(a == L, u, bc)

# Plot

# Add a function for triangulating the solution given by FEniCS
def mesh2triang(mesh):
    import matplotlib.tri as tri
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

# Convert solution into a valid input for plot_trisurf()
object = u.cpp_object()
mesh2 = object.function_space().mesh()
C = object.compute_vertex_values(mesh2)

# Plot
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_trisurf(mesh2triang(mesh2), C, cmap = cm.plasma)
plt.savefig('heateqFEM.png')
plt.show()
```


2D Steady-State Heat Equation Solver
====================================
This README assumes that the reader is familiar with partial differential equations, the finite element method, and the finite difference method. 

## Problem Statement
This project uses both the finite element and finite difference methods to solve the steady-state heat equation, which uses the following equation: $$-k\nabla^2u(x,y)=q(x,y)$$
where $k$ is the thermal conductivity, $u(x,y)$ is the temperature of the material at location $(x,y)$, and $q(x,y)$ is the heat source.

The problem will be solved for a square domain with dimensions $[0,1]\times[0,1]$, and the boundary conditions are $$u(x,y)=1+x^2+\alpha y^2.$$

## Assumptions
The assumptions for solving the heat equation are:
1. The problem has Dirichlet boundary conditions (the ends of the element are held at fixed values).
2. The mesh spacing is constant.
3. The mesh spacing is the same in the $x$ and $y$ directions ($\Delta x=\Delta y$).
4. The thermal conductivity $k$ is constant.
5. The element is homogeneous.

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
The full code for solving the heat equation using the finite element method is:

```python
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Set parameters
k = 0.1
alpha = 3
nx, ny = 100,100
```

```python
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

```

```python
# Solve
u = Function(V)
solve(a == L, u, bc)
```

```python
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

## Set Up Solutions for Solving the Heat Equation with the Finite Difference Method
The full code for solving the heat equation using the finite difference method is:

```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Set up parameters
num_steps = 100
dx = 1.0/num_steps
n = (num_steps+1)*(num_steps+1)
k = 0.1
alpha = 3

# Boundary is u(x,y) = 1 + x^2 + alpha*y^2
def bc(x, y):
    return 1 + x*x + alpha*y*y
```

```python
# Set up matrix A
nx = int(np.sqrt(n))
A = np.zeros((n, n))

# Set up boundary terms
for i in range(nx):
    # Top and bottom boundaries
    A[i,i] = 1.0
    A[i+nx*(nx-1),i+nx*(nx-1)] = 1.0

    # Left and right boundaries
    A[i*nx,i*nx] = 1.0
    A[nx-1+nx*i,nx-1+nx*i] = 1.0

# Set up inner nodes
for i in range(1,nx-1):
    for j in range(1,nx-1):
        A[i+j*nx,i+j*nx] = -4.0
        A[i+j*nx,i+j*nx-1] = 1.0
        A[i+j*nx,i+j*nx+1] = 1.0
        A[i+j*nx,i+j*nx-nx] = 1.0
        A[i+j*nx,i+j*nx+nx] = 1.0
```

```python
# Set up vector b
factor = -dx*dx/k
sourceval = -k*(2+2*alpha)
nx = int(np.sqrt(n))
b = np.zeros((n,1))

tempy1 = 0.0
tempy2 = 1.0

# Set boundary with Dirichlet conditions
for i in range(nx):
    # Bottom row of boundary, where j=0, 0<=i<nx so idx=i+nx*j=i
    tempx1 = i*dx
    b[i] = bc(tempx1, tempy1)
    
    # Top row of boundary, where j=nx-1, 0<=i<nx so idx=i+nx*j=i+nx*(nx-1)
    idx = i+nx*(nx-1)
    b[idx] = bc(tempx1, tempy2)

    # Left most boundary, where 0<=j<nx, i=0 so idx=i+nx*j=nx*j
    idx = i*nx
    b[idx] = bc(tempy1, tempx1)

    # Right most boundary, where 0<=j<nx, i=nx-1 so idx=i+nx*j=nx-1+nx*j
    idx = nx-1+nx*i
    b[idx] = bc(tempy2, tempx1)

# Set source term for every other node
for i in range(1,nx-1):
    tempx1 = i*dx
    for j in range(1,nx-1):
        idx = i+nx*j
        b[idx] = factor*sourceval
```

```python
sol = np.linalg.solve(A,b)
```

```python
# Set up mesh for plotting
nvis = 101
arr = np.zeros((nvis,nvis))

for i in range(nvis):
    for j in range(nvis):
        idx = i+nvis*j
        arr[i,j] = sol[idx]

x = np.linspace(0.0-dx,1.0+dx,101)
y = np.linspace(0.0-dx,1.0+dx,101)

[X,Y] = np.meshgrid(x,y)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_ylim(0,1)
ax.set_xlim(1,0)
ax.view_init(30,30)
ax.plot_surface(X,Y,arr,cmap=cm.plasma)
plt.savefig('heateqFD.png')
```

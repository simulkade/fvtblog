<!-- 
.. title: A simple finite volume toolbox for Matlab
.. slug: 2014-05-20-a-simple-finite-volume-toolbox-for-matlab
.. date: 2014-05-20 22:35:13 UTC+01:00
.. tags: 
.. category: [FVM, Matlab] 
.. link: 
.. description: 
.. type: text
-->

# An introduction to a Matlab FVM (toy) toolbox
For some reasons, I had to solve a few PDE's including single/multi phase flow in porous media, heat transfer in saturated porous media, multi-component mass transfer, and so on. My job never included the development of numerical method. In fact, I was supposed to come up with simple and useful models for the physical system that I was/am studying. Solving those models could be done in PDE solver of my choice. At the end of the day, I chose to write my own codes in Matlab.

## Why [in general]?
Let me count my reasons here. In case you don't like to listen to me bragging about it, you can skip to the middle of this post.

### Black-boxes
These numerical solvers, in particular [COMSOL][1] always was like a black-box to me. I could not understand its error messages, which was frustrating, and even when it generated nice results (which most of the times happens), I was really unable to have a physical feeling for it. I will talk about the physical feeling later, and don't expect a post about love-making.

### Boundary conditions
For some reason, choosing and handling boundary conditions is not explained in books, or papers, or lectures. The same situation applies to the documentation of PDE solvers. Id it really that difficult for you guys to bring us a general Robin boundary condition?

### Geometry
Honestly, one of the most attractive features of PDE solvers is the graphical pre-processor (with all sort of different mesh generation techniques) and CAD modules. For me, with the simple rectangular or cylindrical geometries of my experimental set-up, there was no need of a fancy pre-processor.

### Learning/Teaching
I'll be doing a lot of teaching soon, so I needed to learn numerical methods. What's better than learning by doing? Now, I can share my experiences and my coode with my students.

## Why [in particular]?
I needed a mass conservative scheme (e.g., [finite volume method][2]), which is implemented in an understandable language (yes, I know. C++ is quite beautiful and elegant and understandable even for a kid with the right genes, but I prefer Matlab), with some flexibility for specifying boundary conditions and changing the physics. I actually found a code. It's called [FiPy][3]. But I was already comfortable with Matlab and, don't tell anyone, I couldn't understand [Python][4] and [NumPy][5]. So having FiPy's syntax in mind I decided to write my own code in Matlab. I think it is in good enough shape to be shared with other FVM/Matlab users.

## Really! Why?
For most of our important PDE's in chemical and petroleum engineering (did I mention that I'm a chemical/petroleum engineer?) we have analytical solutions. Also, most of the experiments that we do in the lab can be modeled with a one dimensional PDE. However, so many curious things happen when we go to a two- or three-dimensional domain. I wanted to have something, like FiPy, to make me able to solve a 1D equation, verify it by comparing it to my analytical solution or experimental data, and switch it to 2D and 3D domains without too many modifications. I have it now.

## What do we solve?
We solve this general form of transient convection-diffusion equation:

$$ \alpha\frac{\partial\phi}{\partial t}+\nabla.\left(\mathbf{u}\phi\right)+\nabla.\left(-D\nabla\phi\right)+\beta\phi=\gamma$$ 

with the following general (Robin) boundary condition:

$$a\nabla\phi.\mathbf{n}+b\phi=c.$$

All of the coefficients can be defined explicitly for each control volume or on the surface of a control volume.

## Where to find/How to use 'the code'?
You can download the code from this github repository (click on the `download zip` button): [FVTool][6]
Or alternatively, if you are on linux (and hopefully you are), use the command

```
git clone https://github.com/simulkade/FVTool
```

### How to start
Start [Matlab][7] (or [Octave][8]), go to the `FVTool` folder, and type 

```matlab
FVToolStartUp 
```

You must see a few messages and finally you should see `FiniteVolumeToolbox has started successfully.` in Matlab command prompt.
With the following command, you can see a short document that introduces you to the code:

```matlab
showdemo FVTdemo
```

If you want to jump into it, you can run the following script to solve a diffusion equation with Dirichlet boundary conditions:

```matlab
clc; clear;
L = 50;  % domain length
Nx = 20; % number of cells
m = createMesh3D(Nx,Nx,Nx, L,L,L);
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:)=1; BC.left.c(:)=1; % Dirichlet for the left boundary
BC.right.a(:) = 0; BC.right.b(:)=1; BC.right.c(:)=0; % right boundary
D_val = 1; % value of the diffusion coefficient
D = createCellVariable(m, D_val); % assign the diffusion coefficient to the cells
D_face = harmonicMean(m, D); % calculate harmonic average of the diffusion coef on the cell faces
Mdiff = diffusionTerm(m, D_face); % matrix of coefficients for the diffusion term
[Mbc, RHSbc] = boundaryCondition(m, BC); % matix of coefficients and RHS vector for the BC
M = Mdiff + Mbc; % matrix of cefficients for the PDE
c = solvePDE(m,M, RHSbc); % send M and RHS to the solver
visualizeCells(m, c); % visualize the results
```

You can find more examples [here][9].


  [1]: http://www.comsol.com
  [2]: http://en.wikipedia.org/wiki/Finite_volume_method
  [3]: http://www.ctcms.nist.gov/fipy/
  [4]: https://www.python.org/
  [5]: http://www.numpy.org/
  [6]: https://github.com/simulkade/FVTool
  [7]: http://www.matlab.com
  [8]: http://www.gnu.org/software/octave/
  [9]: https://github.com/simulkade/FVTool/tree/master/Examples/Tutorial

# hBEDA
Numerical solution to the Boltzmann Equation in Diffusion Approximation (BEDA) for spatially homogeneous systems

Numerical solution of the BEDA for an spatially homogeneous system of quarks and gluons. Follow the instructions to run:

1. Compile

`make`

2. Run the code

`./beda args`

Where `args` determine which configuration we choose for the evolution. We can distiguish two main ways to execute the compiled file, as we discuss in the following sections.

## Run with default initial condition

We can easily execute the compiled program for a system initially populated with gluons with a function distribution $f = f_0 \theta (Q-p)$. $Q$ is the characteristic momentum of the system and it is set by default as $Q=1$. In order to run this configuration just do

`./beda type fzero (aS) (np) (dt) (nx) (tMax) (Nf) (Ba) (pMax)`

Where, if not specified, the elements in parenthesis will be setted by default. Only `type` and `fzero` are mandatory. 
- `type` -> Determines which kernels are implemented. The three accepted options are 'e' (only the $2 \leftrightarrow 2$ included), 'i' (only the $1 \leftrightarrow 2$ included) and 'b' (both kernels are implemented).
- `fzero` -> Initial occupancy of the gluon distribution $f_0$.
- `aS` -> Coupling constant of the theory, $\alpha_s$
- `np` -> Number of elements in the momentum grid.
- `dt` -> Time step.
- `nx` -> Number of elements in the $x$ grid, relevant for the calculation of the $1 \leftrightarrow 2$ kernel.
- `tMax` -> Final time of the evolution.
- `Nf` -> Number of quark flavour. If $N_f=0$ the evolution will run without the participation of quarks.
- `Ba` -> Baryon asymetry. Can be 'T' (true) of anything else/default (false).
- `pMax` -> Maximum momentum of the grid. It is relevant to compute the grid element size.

## Run with arbitrary distribution

It is also posible to initialyze the system with an arbitrary configuration saved in a file `filename.dat` (the `.dat` extension is mandatory). Then, we can execute the evolution with

`./beda filename.dat`

We can find a template in [data_init.dat](data_init.dat) with the mandatory format.

- First line must include the name of the files where data will be saved. If files already exist, new data will be appended at the end.
- Second line contain information of 3 relevan quantities. In order: `Initial time`, `Nf` and `nmin`. `nmin` refers to the number of grid elements with a smaller spacing.
- Finally we need at least two columns, but no more than four with the initial distribution. The data in this columns must be, in order: `momentum grid`, `gluon distribution`, `quark distribution` and `antiquark distribution`. Only the first two are mandatory.


## Use the `setxxx` functions to modify the evolution

We can easily modify the evolution previous to compilation time in the [beda.cpp](src/beda.cpp) file using the `setxxx` method of the `KTRun` class after initialization.


## Saving files

Execution generates 5 different files during the evolution where data is saved.
- `filename.fs.dat` -> Saves function distributions and momentum grid for each time step.
- `filename.Cel.dat` -> Saves $2 \leftrightarrow 2$ and momentum grid for each time step.
- `filename.Cinel.dat` -> Saves $1 \leftrightarrow 2$ and momentum grid for each time step.
- `filename.mac.dat` -> Saves different macroscopic quantities at each time step.
- `filename.time.dat` -> Saves initial configuration of the simulation.

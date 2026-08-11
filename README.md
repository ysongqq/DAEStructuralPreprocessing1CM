Structural Preprocessing Method for Nonlinear Differential-Algebraic Equations Using Linear Symbolic Matrices
====

MATLAB/C++ codes for the numerical experiments in the paper
>Structural Preprocessing Method for Nonlinear Differential-Algebraic Equations Using Linear Symbolic Matrices.

It covers the application of the proposed method to the robotic arm DAE and the transistor amplifier DAE; the substitution and augmentation methods used as baselines in the paper
are provided by [DAEPreprocessingToolbox](https://github.com/OptMist-Tokyo/DAEPreprocessingToolbox)
and are not driven from here.

### Requirements

* MATLAB (R2020a was used for the experiments in the paper)
* Symbolic Math Toolbox
* A C++ compiler supported by MATLAB `mex` (GCC 8.1.0 was used)

### Setup

Clone the repository together with the toolbox submodule using `git clone --recursive` command:

If the repository was cloned without `--recursive`:

```
git submodule update --init
```

Then, in MATLAB, build the MEX file and set the search path:

```matlab
cd DAEStructuralPreprocessing1CM
mex -outdir src src/LMMatrixRankMATLAB.cpp
setup
```

`setup` must be run once per MATLAB session. It adds the submodule and the
directories of this repository to the search path, and reloads the MuPAD
package `daepp` shipped with the toolbox.

### Usage

Apply the proposed method to the robotic arm DAE with `K` links (size `3K+2`):

```matlab
>> roboticArm(1)
```

Apply the proposed method to the transistor amplifier DAE:

```matlab
>> circuit()
```

Both functions return the elapsed time in seconds. They also verify that the system
Jacobian of the transformed DAE is nonsingular and raise an error if it is not.

While the method runs, `jacobianToLinearMatrix` prints the total rank of the
coefficient matrices of the linear symbolic matrix approximating the system
Jacobian, together with the number of those matrices. The two values coincide
exactly when every coefficient matrix has rank one, i.e. when the matrix is
already a 1CM-matrix and no rank-one decomposition is needed.


### Directories

```
src/          the proposed method
problems/     the robotic arm DAE with symbolic physical parameters
experiments/  entry points, one per DAE
overlay/      modified copy of one toolbox file
patches/      diff of that file against the upstream version
external/     DAEPreprocessingToolbox (git submodule)
```

### Licence

[MIT License](LICENSE)

### References

> Oki, T., 2023. Improved structural methods for nonlinear differential-algebraic equations via combinatorial relaxation. IMA Journal of Numerical Analysis 43, 357–386.

> De Luca, A., 1988. Control properties of robot arms with joint elasticity.  Analysis and Control of Nonlinear Systems , 61–70.
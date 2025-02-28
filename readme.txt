how to  complie Darwin on personal pc

1. copy mesh_adaptation from megamind, and add export of darwin, prepro, uhmesh in bashrc file
2. download and install openmpi 2.0, add path to .bashrc (add only lib, not bin since bin include mpi, which will conflict with /usr/bin/mpi)
////// Since the original version of darwin is compiled with openmpi 2.0. If you want to use new version openmpi, you have to compile prepro and darwin with the new openmpi.
//////  I do not know how to compile prepro and darwin, maybe need someone who know fortan since the cose is written using fortran.
3. copy 2nd.sh, IN.cfg,spkig.cfg, initial_grid to local folder, modify 2nd.sh slightly (delete slume), and run using "./2nd.sh"

care! 
1. darwin read flow variables as order of d(1), M(6), P(3), T(4), see darwi/code/nerror_estimator.f90
2. but SU2 output as order of d(1), P(3), T(4), M(5), cp(6), see SU2_CFD/src/output/CFlowCompOutput.cpp, line 236-239, line 312-317
3. to make SU2 and darwin consisent, modify SU2 output order, i.e., swap line 238-239, line 314-316, recompile SU2



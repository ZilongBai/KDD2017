# KDD2017 Unsupervised Network Discovery for Brain Imaging Data
This repository contains matlab source code to solve the Orthogonal Nonnegative Matrix tri-Factorization with Spatial Continuity Regularization to learn a complete simplified network with both Nodes and Edges from the corrlation graph of fMRI scan. That is, discover clusters of voxels as Nodes and the inter-/intra-cluster correlations as Edges. Each learnt cluster is spatially continuous. 

The details of the model and the solver are described in our paper
[Unsupervised Network Discovery for Brain Imaging Data](http://dl.acm.org/citation.cfm?id=3098023&CFID=796408940&CFTOKEN=92880021), which demonstrates their utility in the analyses of both synthetic networks constructed with anatomical properties and real-world brain scans. 
The deduction of solving Node discovery subproblem involves techniques similar to Orthogonal Nonnegative Matrix Factorization (Ding et al. 2006). The details of deductions can also be found in our paper.

File(s) in this repository:
- Solver.m: The complete network discovery solver for both the Node discovery subproblem with multiplicative update rules and the Edge discovery problem with nonnegative least squares.
- Theta.m: Calculates Theta matrix according to the spatial coordinates of the vertices of the input correlation graph X.

Please report any bugs/issues to zlbai@ucdavis.edu

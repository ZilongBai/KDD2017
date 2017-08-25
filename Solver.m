function [F,M,ERR,Obj,Conti] = Solver(X,pc,k,MAX,epsilon,Theta);
%
% Solves both the Unsupervised Node Discovery subproblem and the Edge Discovery subproblem (by Zilong Bai, KDD Lab @ University of California, Davis)
%
% The Unsupervised Node Discovery subproblem is solved with multiplicative update rules for F and G. 
%% The detailed formulation and deductions can be found in our paper Unsupervised Network Discovery for Brain Imaging Data, KDD'17. 
%%% The deductions invlolve techniques similar to Orthogonal Nonnegative Matrix Factorization (Ding et al. 2006).
%
% The Edge Discovery subproblem is solved with nonnegative least squares.
%
% Input
%   X: N x N matrix where N is the number of voxels. 
%      Absolute Pearson Correlation graph for all pairs of voxels according to their temporal patterns. 
%      In the paper we evaluated all voxels within the brain region on one slice (36th).
%   pc: Scalar. Weight of the continuity regularization term.
%   k: Integer. Number of Nodes to be discover
%   MAX: Integer.  maximum number of iterations
%   epsilon: Scalar. An extremely small postive scalar to replace non-positive elements to avoid multiplicative update rules being trapped by zero values.
%   Theta: N x N matrix where N is the number of voxels. Penalty matrix calculated according to the spatial coordinates of the voxels.
%
% Output
%   F: N x k matrix.
%      Cluster indicator matrix for discovered Nodes. Each Node is inidcated by a column vector of F.
%   M: k x k matrix.
%      Inter-/intra-cluster correlation matrix for discovered Edges.
%   ERR, Obj, Conti: Vectors for debugging purpose. 
%                    Record the reconstruction error, objective function value, and the continuity regularization term value of each iteration while solving the Unsupervised Node discovery subproblem

[spatio, spatioy] = size(X); % X is a symmetric absolute Pearson correlation matrix. Hence spatio equals spatioy, denoting the number of vertices on this correlation graph, i.e. the voxels in the spatial domain. Note they are NOT the spatial coordinates of the voxels.

%Gaussian Kernel Theta for spatial continuity regularization can be pre-calculated.

%Initialization

F = rand(spatio, k);

G = rand(spatioy, k);

%=== Node Discovery Subproblem

% Start alternative solving process for F and G. F is the indicator matrix for Nodes.

ERR = [];

for iteration = 1:MAX

% Updating F with multiplicative update rule

SQF = ((X*G)./(F*(G'*G)+pc*Theta*F + F*(F'*X*G - (G'*G) -F'*pc*Theta*F))); % Lambda in the paper is substituted with its explicit expression with F, G, pc, and Theta. 
SQF(find(SQF<=0)) = epsilon; 
F = F.*(sqrt(SQF));
%F(find(F<=0)) = epsilon;

% Updating G with multiplicative update rule
G = G.*((X'*F)./(G*(F'*F)));
G(find(G<=0)) = epsilon;

% Normalizing columns of F
for col = 1:k
	normF = norm(F(:,col));
        F(:,col) = F(:,col)./normF;
	G(:,col) = G(:,col).*normF;
end

ERR(iteration) = norm(X-F*G','fro')/norm(X,'fro');
Obj(iteration) = (norm(X-F*G','fro'))^2 + trace((F')*Theta*F);
Conti(iteration) = trace(F'*Theta*F); 
end

%=== Edge Discovery Subproblem

% Solving symmetrical M with Kronecker product and Nonnegative Least Squares. M records the Edges. 

vecX = reshape(X,spatio*spatioy,1);

kronF = kron(F,F);

vecM = lsqnonneg(kronF,vecX);

M = reshape(vecM,k,k);


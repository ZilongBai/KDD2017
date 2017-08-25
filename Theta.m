function [THETA] = Theta(X,Y,sigma)
%
% Calculates the spatial continuity penalty matrix THETA according to the spatial coordinates of the voxels. The penalty for each pair of voxels being clustered into one cluster is in the form of reciprocal Gaussian Kernel. The mathematical formulation can be founed in our paper Unsupervised Network Discovery for Brain Imaging Data, KDD'17
% (by Zilong Bai, KDD Lab @ University of California, Davis)
%
% Input:
%    X: N x s matrix, where N is the number of voxels while s in the dimension of the coordinates of each voxel in the spatial domain.
%    Y: N x s matrix, where N is the number of voxels while s in the dimension of the coordinates of each voxel in the spatial domain.
%    In this paper X equals Y. This function can be used for more general settings.
%    sigma: bandwidth for the penalty kernel calculation. It impacts the size of the continuous Nodes to be discovered. 
% Output:
%    THETA: N x N penalty matrix for each pair of voxels.

[xx,xy] = size(X);
[yx,yy] = size(Y);

XY = [];
for i = 1:xx
for j = 1:yx
  XY(i,j) = (norm(X(i,:)-Y(j,:)))^2;
end
end

THETA = exp(XY./(2*sigma^2)); 

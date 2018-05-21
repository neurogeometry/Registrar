% This function finds optimal linear (affine) transformation (X -> X') for 
% assignment Y=L*X+b. X and Y are 3xN
% Transformation is regularized with energy density of elastic deformation 

function [L,b]=Optimal_Affine_Transform(X,Y)

% regularization parameter
mu=10^4;

N=size(X,2);
X_cm=X-mean(X,2)*ones(1,N);
Y_cm=Y-mean(Y,2)*ones(1,N);
Cov_XX=X_cm*X_cm'./N;
Cov_YX=Y_cm*X_cm'./N;

L=(Cov_YX+mu.*diag(ones(1,size(X,1))))/(Cov_XX+mu.*diag(ones(1,size(X,1))));
b=mean(Y,2)-mean(L*X,2);
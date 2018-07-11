% This function finds optimal linear (affine) transformation (X -> X') for 
% assignment Y=L*X+b. X and Y are 3xN
% Transformation is regularized with energy density of elastic deformation 

function [X_aligned,L,b]=Optimal_Affine_Transform(X,Y,mu)

% regularization parameter
% mu=0.1;

N=size(X,2);
X_cm=X-mean(X,2)*ones(1,N);
Y_cm=Y-mean(Y,2)*ones(1,N);
Cov_XX=X_cm*X_cm'./N;
Cov_YX=Y_cm*X_cm'./N;

L=(Cov_YX+mu.*diag(ones(1,size(X,1))))/(Cov_XX+mu.*diag(ones(1,size(X,1))));
b=mean(Y,2)-mean(L*X,2);

X_aligned=L*X+b*ones(1,size(X,2));
% disp(([mean(sum((Y-X).^2)), mean(sum((Y-X_aligned).^2))]).^0.5)
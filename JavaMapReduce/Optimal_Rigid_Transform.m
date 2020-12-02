% This function finds optimal rigid (translation+rotation) transformation (X -> X') for 
% assignment Y=R*X+b. X and Y are 3xN
% Solution is based on svd of covariance matrix 

function [R,b]=Optimal_Rigid_Transform(X,Y)

N=size(X,2);
X_cm=X-mean(X,2)*ones(1,N);
Y_cm=Y-mean(Y,2)*ones(1,N);
Cov_XY=X_cm*Y_cm'./N;
[U,~,V] = svd(Cov_XY);

R=V*U';
b=mean(Y,2)-mean(R*X,2);

% X_aligned=R*X+b*ones(1,size(X,2));
% disp(([mean(sum((Y-X).^2)), mean(sum((Y-X_aligned).^2))]).^0.5)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(1,2,1)
% plot3(X(1,:),X(2,:),X(3,:),'r*')
% axis equal, hold on, box on
% plot3(Y(1,:),Y(2,:),Y(3,:),'k*')
% title('Original')
% subplot(1,2,2)
% plot3(X_aligned(1,:),X_aligned(2,:),X_aligned(3,:),'r*')
% axis equal, hold on, box on
% plot3(Y(1,:),Y(2,:),Y(3,:),'k*')
% title('Aligned')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % visualize transformation
% addpath('C:\Armen\DIADEM\Neuron Tracer V20')
% nx=10; ny=10; nz=10;
% Min=min([min(X,[],2),min(Y,[],2)],[],2);
% Max=max([max(X,[],2),max(Y,[],2)],[],2);
% [xx,yy,zz]=ndgrid(1:nx,1:ny,1:nz);
% r_grid=[Min(1)+(Max(1)-Min(1)).*(xx(:)-1)./(nx-1),Min(2)+(Max(2)-Min(2)).*(yy(:)-1)./(ny-1),Min(3)+(Max(3)-Min(3)).*(zz(:)-1)./(nz-1)]';
% AM_grid=(abs(bsxfun(@minus,xx(:),xx(:)'))+abs(bsxfun(@minus,yy(:),yy(:)'))+abs(bsxfun(@minus,zz(:),zz(:)')))==1;
% figure
% subplot(1,2,1)
% PlotAM(AM_grid,r_grid')
% title('Original')
% N=nx*ny*nz;
% r_grid_aligned=R*r_grid+b*ones(1,N);
% subplot(1,2,2)
% PlotAM(AM_grid,r_grid_aligned')
% title('Aligned')
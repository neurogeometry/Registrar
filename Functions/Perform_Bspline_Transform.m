% This function performs a Cubic B-spline transform X'=X+f(X)
% X and X' are 3xN. E is the deformation energy
% Nxyz is a vector defining the number of a B-spline cells
% nxyz is a vector defining the size of a B-spline cell
% Cxyz is a 3x(prod(Nxyz + 3)) matrics of B-spline coefficients

function [X_prime,E]=Perform_Bspline_Transform(X,Cxyz,Nxyz,nxyz)

NB=4; % # of B-spline polynomials

[ll,mm,nn]=ind2sub([NB,NB,NB],(1:NB^3));
[ii,jj,kk]=ind2sub(Nxyz,1:prod(Nxyz));

B=@(x) cat(1,(1-x).^3./6,(3.*x.^3-6.*x.^2+4)./6,(-3.*x.^3+3.*x.^2+3.^x+1)./6,x.^3./6);

C=[1/252,43/1680,1/84,1/5040;
   43/1680,33/140,311/1680,1/84;
   1/84,311/1680,33/140,43/1680;
   1/5040,1/84,43/1680,1/252];
D=[1/20,7/120,-1/10,-1/120;
   7/120,17/60,-29/120,-1/10;
   -1/10,-29/120,17/60,7/120;
   -1/120,-1/10,7/120,1/20];
Q=(D(ll,ll).*C(mm,mm).*C(nn,nn)./nxyz(1)^2+C(ll,ll).*D(mm,mm).*C(nn,nn)./nxyz(2)^2+C(ll,ll).*C(mm,mm).*D(nn,nn)./nxyz(3)^2)./prod(Nxyz);

X_prime=X;
E=0;
for i=1:length(ii)
    ind=(ii(i)-1<=X(1,:)./nxyz(1) & X(1,:)./nxyz(1)<=ii(i) & jj(i)-1<=X(2,:)./nxyz(2) & X(2,:)./nxyz(2)<=jj(i) & kk(i)-1<=X(3,:)./nxyz(3) & X(3,:)./nxyz(3)<=kk(i));
    Bx=B(X(1,ind)./nxyz(1)-ii(i)+1);
    By=B(X(2,ind)./nxyz(2)-jj(i)+1);
    Bz=B(X(3,ind)./nxyz(3)-kk(i)+1);
    temp_ind=sub2ind(Nxyz+NB-1,ii(i)-1+ll,jj(i)-1+mm,kk(i)-1+nn);
    CC=Cxyz(:,temp_ind);
    
    X_prime(:,ind)=X_prime(:,ind)+CC*(Bx(ll,:).*By(mm,:).*Bz(nn,:));
    E=E+CC(1,:)*Q*CC(1,:)'+CC(2,:)*Q*CC(2,:)'+CC(3,:)*Q*CC(3,:)';
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % visualize transformation
% addpath('C:\Armen\DIADEM\Neuron Tracer V20')
% nx=11; ny=11; nz=11;
% [xx,yy,zz]=ndgrid(1:nx,1:ny,1:nz);
% X=[(xx(:)-1)./(nx-1),(yy(:)-1)./(ny-1),(zz(:)-1)./(nz-1)];
% AM_grid=(abs(bsxfun(@minus,xx(:),xx(:)'))+abs(bsxfun(@minus,yy(:),yy(:)'))+abs(bsxfun(@minus,zz(:),zz(:)')))==1;
% figure
% subplot(1,2,1)
% PlotAM(AM_grid,X)
% title('Original')
% 
% subplot(1,2,2)
% PlotAM(AM_grid,X_prime)
% title('Aligned')
% 
% del_r2=sum((r_grid_aligned-r_grid).^2,1);
% disp([mean(del_r2)^0.5, max(del_r2)^0.5])

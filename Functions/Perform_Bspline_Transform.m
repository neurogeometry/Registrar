% This function performs a Cubic B-spline transform X'=X+f(X)
% X and X' are 3xN. E is the deformation energy
% Nxyz is a vector defining the number of a B-spline cells
% nxyz is a vector defining the size of a B-spline cell
% Cxyz is a 3x(prod(Nxyz + 3)) matrics of B-spline coefficients

function [X_prime,StackPosition_prime,E]=Perform_Bspline_Transform(X,StackPosition,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine)

X_prime=[];
StackPosition_prime=[];
E=[];
sizeX=size(X);

if length(sizeX)==2
    sizeX(3)=1;
end

if isempty(StackPosition)  % X is a 3xN array of positions
    [X_prime,E]=Perform_Bspline_Transform0(X,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine,1);
    
elseif length(sizeX)==3 % X is a 2d image or a 3d image stack
    
    Verts=[1,1,1;1,1,sizeX(3);1,sizeX(2),1;sizeX(1),1,1;1,sizeX(2),sizeX(3);sizeX(1),1,sizeX(3);sizeX(1),sizeX(2),1;sizeX(1),sizeX(2),sizeX(3)]';
    [Verts_prime,~]=Perform_Bspline_Transform0(Verts,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine,0); 
    Min=round(min(Verts_prime,[],2));
    Max=round(max(Verts_prime,[],2));
    sizeX_prime=[Max(1)-Min(1)+1,Max(2)-Min(2)+1,Max(3)-Min(3)+1];
    X_prime=zeros(sizeX_prime,class(X));
    
    [xx,yy,zz]=ndgrid(1:sizeX(1),1:sizeX(2),1:sizeX(3));
    xyz=[xx(:),yy(:),zz(:)]';
    clear xx yy zz
    [xyz_prime,~]=Perform_Bspline_Transform0(xyz,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine,0);
    xyz_prime=round(xyz_prime);
    xyz_prime=xyz_prime-(Min-1)*ones(1,prod(sizeX));
    
    ind=(xyz_prime(1,:)>=1 & xyz_prime(1,:)<=sizeX_prime(1) & xyz_prime(2,:)>=1 & xyz_prime(2,:)<=sizeX_prime(2) & xyz_prime(3,:)>=1 & xyz_prime(3,:)<=sizeX_prime(3));
    ind2=xyz_prime(1,ind)+(xyz_prime(2,ind)-1).*sizeX_prime(1)+(xyz_prime(3,ind)-1).*(sizeX_prime(1)*sizeX_prime(2));
    X_prime(ind2)=X(ind);
    
    StackPosition_prime=StackPosition+Min-1;
      
    %     figure
    %     imshow(max(X,[],3),[0 max(X(:))])
    %     figure
    %     imshow(max(X_prime,[],3),[0 max(X(:))])
    
%     xyz_prime0=xyz;
%     xyz_prime=xyz_prime0;
%     for i=1:10
%         [xyz_prime_temp,~]=Perform_Bspline_Transform0(xyz_prime,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine,0);
%         xyz_prime=xyz_prime-0.01.*(xyz_prime_temp-xyz_prime0);
%         (mean(sum((xyz_prime-xyz_prime0).^2))).^0.5
%     end
    
    
else
    disp('Format of X is incorrect.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_prime,E]=Perform_Bspline_Transform0(X,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine,energy)

NB=4; % # of B-spline polynomials

[ll,mm,nn]=ind2sub([NB,NB,NB],(1:NB^3));
[ii,jj,kk]=ind2sub(Nxyz,1:prod(Nxyz));

B=@(x) cat(1,(1-x).^3./6,(3.*x.^3-6.*x.^2+4)./6,(-3.*x.^3+3.*x.^2+3.*x+1)./6,x.^3./6);

if affine==1
    [X,~]=Perform_Linear_Transform(X',[],L,b);
    X=X';
end
X=X-Grid_start*ones(1,size(X,2));

if energy==1
    C=[1/252,43/1680,1/84,1/5040;
        43/1680,33/140,311/1680,1/84;
        1/84,311/1680,33/140,43/1680;
        1/5040,1/84,43/1680,1/252];
    D=[1/20,7/120,-1/10,-1/120;
        7/120,17/60,-29/120,-1/10;
        -1/10,-29/120,17/60,7/120;
        -1/120,-1/10,7/120,1/20];
    Q=(D(ll,ll).*C(mm,mm).*C(nn,nn)./nxyz(1)^2+C(ll,ll).*D(mm,mm).*C(nn,nn)./nxyz(2)^2+C(ll,ll).*C(mm,mm).*D(nn,nn)./nxyz(3)^2)./prod(Nxyz);
    E=0;
else
   E=[]; 
end

X_prime=X;
for i=1:length(ii)
    ind=((ii(i)-1<=X(1,:)./nxyz(1)) & (X(1,:)./nxyz(1)<ii(i)) & (jj(i)-1<=X(2,:)./nxyz(2)) & (X(2,:)./nxyz(2)<jj(i)) & (kk(i)-1<=X(3,:)./nxyz(3)) & (X(3,:)./nxyz(3)<kk(i)));
    Bx=B(X(1,ind)./nxyz(1)-ii(i)+1);
    By=B(X(2,ind)./nxyz(2)-jj(i)+1);
    Bz=B(X(3,ind)./nxyz(3)-kk(i)+1);
    temp_ind=sub2ind(Nxyz+NB-1,ii(i)-1+ll,jj(i)-1+mm,kk(i)-1+nn);
    CC=Cxyz(:,temp_ind);
    
    X_prime(:,ind)=X_prime(:,ind)+CC*(Bx(ll,:).*By(mm,:).*Bz(nn,:));
    if energy==1
        E=E+CC(1,:)*Q*CC(1,:)'+CC(2,:)*Q*CC(2,:)'+CC(3,:)*Q*CC(3,:)';
    end
end
X_prime=X_prime+Grid_start*ones(1,size(X,2));

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

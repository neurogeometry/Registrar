% This function finds optimal non-rigid transformation (X -> X') based on 
% Legendre expansion for assignment Y=f(X), where X and Y are 3xN.
% X' = X + sum(X_{lmn} P_l P_m P_n)
% N_L is the number of Legendre pilynomials used
% The problem is regularized with elastic deformation energy, assuming 
% linear, homogeneous, and isotropic elasticity

function [XYZlmn,X_aligned]=Optimal_Nonrigid_Transform_org(X,Y,N_L)

lambda0=0.0001; % regularizer

% Lame parameters
% lambda=0;
% mu=0.1;

addpath('C:\Armen\DIADEM\Neuron Tracer V20')
N=size(X,2);

if N_L==0
    P = @(x) cat(3,zeros(size(x)));
elseif N_L==1
    P = @(x) cat(3,ones(size(x)));
elseif N_L==2
    P = @(x) cat(3,ones(size(x)),x);
elseif N_L==3
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2);
elseif N_L==4
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2,(5.*x.^3-3.*x)./2);
elseif N_L==5
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2,(5.*x.^3-3.*x)./2,(35.*x.^4-30.*x.^2+3)./8);
elseif N_L==6
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2,(5.*x.^3-3.*x)./2,(35.*x.^4-30.*x.^2+3)./8,...
     (63.*x.^5-70.*x.^3+15.*x)./8);
elseif N_L==7
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2,(5.*x.^3-3.*x)./2,(35.*x.^4-30.*x.^2+3)./8,...
        (63.*x.^5-70.*x.^3+15.*x)./8,(231.*x.^6-315.*x.^4+105.*x.^2-5)./16);
elseif N_L==8
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2,(5.*x.^3-3.*x)./2,(35.*x.^4-30.*x.^2+3)./8,...
        (63.*x.^5-70.*x.^3+15.*x)./8,(231.*x.^6-315.*x.^4+105.*x.^2-5)./16,...
        (6435.*x.^8-12012.*x.^6+6930.*x.^4-1260.*x.^2+35)./128);
elseif N_L==9
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2,(5.*x.^3-3.*x)./2,(35.*x.^4-30.*x.^2+3)./8,...
        (63.*x.^5-70.*x.^3+15.*x)./8,(231.*x.^6-315.*x.^4+105.*x.^2-5)./16,...
        (6435.*x.^8-12012.*x.^6+6930.*x.^4-1260.*x.^2+35)./128,...
        (12115.*x.^9-25740.*x.^7+18018.*x.^5-4620.*x.^3+315.*x)./128);
elseif N_L==10
    P = @(x) cat(3,ones(size(x)),x,(3.*x.^2-1)./2,(5.*x.^3-3.*x)./2,(35.*x.^4-30.*x.^2+3)./8,...
        (63.*x.^5-70.*x.^3+15.*x)./8,(231.*x.^6-315.*x.^4+105.*x.^2-5)./16,...
        (6435.*x.^8-12012.*x.^6+6930.*x.^4-1260.*x.^2+35)./128,...
        (12115.*x.^9-25740.*x.^7+18018.*x.^5-4620.*x.^3+315.*x)./128,...
        (46189.*x.^10-109395.*x.^8+90090.*x.^6-30030.*x.^4+3465.*x.^2-63)./256);
end

% squeeze the traces into -1:1 range
Min=min([min(X,[],2),min(Y,[],2)],[],2);
Max=max([max(X,[],2),max(Y,[],2)],[],2);
X_scaled=(X-Min*ones(1,N))./(((Max-Min)./2*ones(1,N)))-1;
Y_scaled=(Y-Min*ones(1,N))./(((Max-Min)./2*ones(1,N)))-1;

[l,m,n]=ind2sub([N_L,N_L,N_L],1:N_L^3);
Plmn=zeros(N,N_L^3);
PP=P(Y_scaled);
for i=1:N_L^3
    Plmn(:,i)=PP(1,:,l(i)).* PP(2,:,m(i)).*PP(3,:,n(i));
    %Plmn(:,i)=legendreP(l(i)-1,Y_scaled(1,:)).*legendreP(m(i)-1,Y_scaled(2,:)).*legendreP(n(i)-1,Y_scaled(3,:));
end


% with simple regularization
XYZlmn=((1/N).*Plmn'*Plmn+(lambda0*N_L^3).*diag(ones(1,N_L^3)))\(-(N_L^3/N).*(X_scaled-Y_scaled)*Plmn)';

%{
% with reqularization based on elastic deformation energy
L=(Max-Min)./2;
XX=sparse(N_L^3,N_L^3);
YY=sparse(N_L^3,N_L^3);
ZZ=sparse(N_L^3,N_L^3);
YX=sparse(N_L^3,N_L^3);
XY=sparse(N_L^3,N_L^3);
ZX=sparse(N_L^3,N_L^3);
XZ=sparse(N_L^3,N_L^3);
YZ=sparse(N_L^3,N_L^3);
ZY=sparse(N_L^3,N_L^3);
for i=1:N_L^3
    for j=1:N_L^3
        XX(i,j)=2*min([l(i),l(j)])*(min([l(i),l(j)])+1)*(1+(-1)^(l(i)+l(j)))*(m(i)==m(j))/(2*m(i)+1)*(n(i)==n(j))/(2*n(i)+1);
        YY(i,j)=2*min([m(i),m(j)])*(min([m(i),m(j)])+1)*(1+(-1)^(m(i)+m(j)))*(l(i)==l(j))/(2*l(i)+1)*(n(i)==n(j))/(2*n(i)+1);
        ZZ(i,j)=2*min([n(i),n(j)])*(min([n(i),n(j)])+1)*(1+(-1)^(n(i)+n(j)))*(l(i)==l(j))/(2*l(i)+1)*(m(i)==m(j))/(2*m(i)+1);
        YX(i,j)=2*(l(j)>l(i))*(1-(-1)^(l(i)+l(j)))*(m(i)>m(j))*(1-(-1)^(m(i)+m(j)))*(n(i)==n(j))/(2*n(i)+1);
        XY(i,j)=2*(l(i)>l(j))*(1-(-1)^(l(i)+l(j)))*(m(j)>m(i))*(1-(-1)^(m(i)+m(j)))*(n(i)==n(j))/(2*n(i)+1);
        ZX(i,j)=2*(l(j)>l(i))*(1-(-1)^(l(i)+l(j)))*(n(i)>n(j))*(1-(-1)^(n(i)+n(j)))*(m(i)==m(j))/(2*m(i)+1);
        XZ(i,j)=2*(l(i)>l(j))*(1-(-1)^(l(i)+l(j)))*(n(j)>n(i))*(1-(-1)^(n(i)+n(j)))*(m(i)==m(j))/(2*m(i)+1);
        YZ(i,j)=2*(n(j)>n(i))*(1-(-1)^(n(i)+n(j)))*(m(i)>m(j))*(1-(-1)^(m(i)+m(j)))*(l(i)==l(j))/(2*l(i)+1);
        ZY(i,j)=2*(n(i)>n(j))*(1-(-1)^(n(i)+n(j)))*(m(j)>m(i))*(1-(-1)^(m(i)+m(j)))*(l(i)==l(j))/(2*l(i)+1);
    end
end
B=(1/N).*Plmn'*Plmn;
A=(1/16).*[16.*B+((2*mu+lambda)./L(1)^2).*XX+(mu/L(2)^2).*YY+(mu/L(3)^2).*ZZ,(mu/L(1)/L(2)).*YX+(lambda/L(1)/L(2)).*XY,(mu/L(1)/L(3)).*ZX+(lambda/L(1)/L(3)).*XZ;...
    (mu/L(1)/L(2)).*XY+(lambda/L(1)/L(2)).*YX,16.*B+(mu/L(1)^2).*XX+((2*mu+lambda)/L(2)^2).*YY+(mu/L(3)^2).*ZZ,(mu/L(2)/L(3)).*ZY+(lambda/L(2)/L(3)).*YZ;...
    (mu/L(1)/L(3)).*XZ+(lambda/L(1)/L(3)).*ZX,(mu/L(2)/L(3)).*YZ+(lambda/L(2)/L(3)).*ZY,16.*B+(mu/L(1)^2).*XX+(mu/L(2)^2).*YY+((2*mu+lambda)/L(3)^2).*ZZ];
b=(-(N_L^3/N).*(X_scaled-Y_scaled)*Plmn)';
XYZlmn=zeros(size(b));
XYZlmn(:)=A\b(:);
%}

X_aligned=X_scaled+(1/N_L^3).*XYZlmn'*Plmn'; 
X_aligned=(X_aligned+1).*(((Max-Min)./2*ones(1,N)))+Min*ones(1,N);
disp(([mean(sum((Y-X).^2)), mean(sum((Y-X_aligned).^2))]).^0.5)


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
% nx=10; ny=10; nz=10;
% [xx,yy,zz]=ndgrid(1:nx,1:ny,1:nz);
% r_grid=[Min(1)+(Max(1)-Min(1)).*(xx(:)-1)./(nx-1),Min(2)+(Max(2)-Min(2)).*(yy(:)-1)./(ny-1),Min(3)+(Max(3)-Min(3)).*(zz(:)-1)./(nz-1)]';
% AM_grid=(abs(bsxfun(@minus,xx(:),xx(:)'))+abs(bsxfun(@minus,yy(:),yy(:)'))+abs(bsxfun(@minus,zz(:),zz(:)')))==1;
% figure
% subplot(1,2,1)
% PlotAM(AM_grid,r_grid')
% title('Original')
% 
% N=nx*ny*nz;
% r_grid_scaled=(r_grid-Min*ones(1,N))./(((Max-Min)./2*ones(1,N)))-1;
% Plmn=zeros(N,N_L^3);
% PP=P(r_grid_scaled);
% for i=1:N_L^3
%     Plmn(:,i)=PP(1,:,l(i)).* PP(2,:,m(i)).*PP(3,:,n(i));
% end
% r_grid_aligned=r_grid_scaled+(1/N_L^3).*XYZlmn'*Plmn';
% r_grid_aligned=(r_grid_aligned+1).*(((Max-Min)./2*ones(1,N)))+Min*ones(1,N);
% subplot(1,2,2)
% PlotAM(AM_grid,r_grid_aligned')
% title('Aligned')
% 
% del_r2=sum((r_grid_aligned-r_grid).^2,1);
% disp([mean(del_r2)^0.5, max(del_r2)^0.5])

% This function finds optimal Cubic affine + B-spline transform X'=LX+b, Y=X'+f(X')
% X and Y are 3xN
% Transformation is regularized with energy density of elastic deformation 
% L, b, Cxyz, Nxyz, nxyz, and Grid_start defines the transformation
% affine = 1/0 turns affine transform on/off

function [X_aligned,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(X,Y,nxyz,affine,mu)

if affine==1
    [X_affine,L,b]=Optimal_Affine_Transform(X,Y,mu);
elseif affine==0
    X_affine=X;
    L=[];
    b=[];
end
% if isnan(X_affine) 
%     X_affine=X;
%     L=[];
%     b=[];
%     affine=0;
% end

Min=min(X_affine,[],2);
Max=max(X_affine,[],2);
Nxyz=ceil((Max-Min+1)./nxyz);
Grid_start=(Min+Max-Nxyz.*nxyz)./2;

X_affine=X_affine-Grid_start*ones(1,size(X,2));
Y=Y-Grid_start*ones(1,size(Y,2));

% regularization parameter
% mu=1; %[0, 2.^[-10:10]]
learning_rate=5;
Max_iterations=100000;
TolCost=10^-6;

N=size(X,2);
NB=4; % # of B-spline polynomials

[ll,mm,nn]=ind2sub([NB,NB,NB],(1:NB^3)');
[il,jm,kn]=ind2sub(Nxyz+NB-1,1:prod(Nxyz+NB-1));

B=@(x) cat(1,(1-x).^3./6,(3.*x.^3-6.*x.^2+4)./6,(-3.*x.^3+3.*x.^2+3.*x+1)./6,x.^3./6);

C=[1/252,43/1680,1/84,1/5040;
   43/1680,33/140,311/1680,1/84;
   1/84,311/1680,33/140,43/1680;
   1/5040,1/84,43/1680,1/252];
D=[1/20,7/120,-1/10,-1/120;
   7/120,17/60,-29/120,-1/10;
   -1/10,-29/120,17/60,7/120;
   -1/120,-1/10,7/120,1/20];
Q=(D(ll,ll).*C(mm,mm).*C(nn,nn)./nxyz(1)^2+C(ll,ll).*D(mm,mm).*C(nn,nn)./nxyz(2)^2+C(ll,ll).*C(mm,mm).*D(nn,nn)./nxyz(3)^2)*2./prod(Nxyz);

dfdXxyz=zeros(length(il),N);
for i=1:length(il)
    for j=1:NB^3
        if 0<=il(i)-ll(j) && il(i)-ll(j)<=Nxyz(1)-1 && 0<=jm(i)-mm(j) && jm(i)-mm(j)<=Nxyz(2)-1 && 0<=kn(i)-nn(j) && kn(i)-nn(j)<=Nxyz(3)-1
            tempx=(il(i)-ll(j)<=X_affine(1,:)./nxyz(1) & X_affine(1,:)./nxyz(1)<=il(i)-ll(j)+1);
            tempy=(jm(i)-mm(j)<=X_affine(2,:)./nxyz(2) & X_affine(2,:)./nxyz(2)<=jm(i)-mm(j)+1);
            tempz=(kn(i)-nn(j)<=X_affine(3,:)./nxyz(3) & X_affine(3,:)./nxyz(3)<=kn(i)-nn(j)+1);
            ind=(tempx & tempy & tempz);
            if any(ind)
                Bx=B(X_affine(1,ind)./nxyz(1)-il(i)+ll(j));
                By=B(X_affine(2,ind)./nxyz(2)-jm(i)+mm(j));
                Bz=B(X_affine(3,ind)./nxyz(3)-kn(i)+nn(j));
                dfdXxyz(i,ind)=dfdXxyz(i,ind)+Bx(ll(j),:).*By(mm(j),:).*Bz(nn(j),:);
            end
        end
    end
end

count=0;
E=0;
Cost=mean(sum((Y-X_affine).^2));
delCost=inf;
Cxyz=zeros(3,prod(Nxyz + 3));
X_aligned=X_affine;
while count<=Max_iterations && (delCost<-TolCost || delCost>0)
    count=count+1;
%      disp([E, mean(sum((Y-X_aligned).^2)).^0.5, Cost])
    
    dEdXYZ=zeros(3,prod(Nxyz+NB-1));
    for i=1:length(il)
        S=0;
        for j=1:NB^3
            if 0<=il(i)-ll(j) && il(i)-ll(j)<=Nxyz(1)-1 && 0<=jm(i)-mm(j) && jm(i)-mm(j)<=Nxyz(2)-1 && 0<=kn(i)-nn(j) && kn(i)-nn(j)<=Nxyz(3)-1
                temp_ind=(il(i)-ll(j)+ll)+(jm(i)-mm(j)+mm-1).*(Nxyz(1)+NB-1)+(kn(i)-nn(j)+nn-1).*(Nxyz(1)+NB-1)*(Nxyz(2)+NB-1);
                S=S+Cxyz(:,temp_ind)*Q(:,j);
            end
        end
        dEdXYZ(:,i)=dEdXYZ(:,i)+(2/prod(Nxyz)).*S;
    end
    
    Cxyz_new=Cxyz-learning_rate.*((2/N).*(X_aligned-Y)*dfdXxyz'+mu.*dEdXYZ);
    [X_aligned_new,~,E_new]=Perform_Bspline_Transform(X_affine,[],[],[],Cxyz_new,Nxyz,nxyz,[0;0;0],0); 
    Cost_new=mean(sum((Y-X_aligned_new).^2))+mu*E_new;
    
    delCost=2*(Cost_new-Cost)/(Cost_new+Cost);
    if delCost<0
        X_aligned=X_aligned_new;
        E=E_new;
        Cxyz=Cxyz_new;
        Cost=Cost_new;    
    else
        learning_rate=learning_rate/1.1;
%          disp(learning_rate)
    end
end
% count
% learning_rate
X_aligned=X_aligned+Grid_start*ones(1,size(X,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y=Y+Grid_start*ones(1,size(Y,2));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualize transformation
% % addpath('C:\Armen\DIADEM\Neuron Tracer V20')
% nx=10; ny=10; nz=10;
% Min=min(X,[],2);
% Max=max(X,[],2);
% [xx,yy,zz]=ndgrid(1:nx,1:ny,1:nz);
% r_grid=[Min(1)+(Max(1)-Min(1)).*(xx(:)-1)./(nx-1),Min(2)+(Max(2)-Min(2)).*(yy(:)-1)./(ny-1),Min(3)+(Max(3)-Min(3)).*(zz(:)-1)./(nz-1)]';
% AM_grid=(abs(bsxfun(@minus,xx(:),xx(:)'))+abs(bsxfun(@minus,yy(:),yy(:)'))+abs(bsxfun(@minus,zz(:),zz(:)')))==1;
% figure
% subplot(1,2,1)
% PlotAM(AM_grid,r_grid','k')
% title('Original')
% [r_grid_aligned,~,~]=Perform_Bspline_Transform(r_grid,[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
% subplot(1,2,2)
% PlotAM(AM_grid,r_grid_aligned','k')
% title('Aligned')

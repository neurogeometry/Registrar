% This function finds optimal non-rigid transformation (X -> X') based on 
% Legendre expansion for assignment Y=f(X), where X and Y are 3xN.
% X' = X + sum(X_{lmn} P_l P_m P_n)
% N_L is the number of Legendre polynomials used
% Min and Max are 3x1 vectors specifying the transformation region. All
% points must lie within the box specified by these vectors.
% The problem is regularized with elastic deformation energy
% For best results, perform affine transform prior to non-rigid

function [X_aligned,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(X,Y,N_L,Min,Max)


% Regularization parameter
mu=10^4;

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

N=size(X,2);
X_scaled=(2.*X-(Min+Max)*ones(1,N))./((Max-Min)*ones(1,N));
Y_scaled=(2.*Y-(Min+Max)*ones(1,N))./((Max-Min)*ones(1,N));

[l,m,n]=ind2sub([N_L,N_L,N_L],1:N_L^3);
Plmn=zeros(N,N_L^3);
PP=P(X_scaled);
for i=1:N_L^3
    Plmn(:,i)=PP(1,:,l(i)).* PP(2,:,m(i)).*PP(3,:,n(i));
end

L=(Max-Min)./2;
XX=zeros(N_L^3,N_L^3);
YY=zeros(N_L^3,N_L^3);
ZZ=zeros(N_L^3,N_L^3);
b=zeros(N_L^3,3);
for i=1:N_L^3
    ind_mn=(m==m(i) & n==n(i));
    ind_ln=(l==l(i) & n==n(i));
    ind_lm=(l==l(i) & m==m(i));
    
    min_l=min(l(i),l(ind_mn));
    min_m=min(m(i),m(ind_ln));
    min_n=min(n(i),n(ind_lm));
    
    XX(i,ind_mn)=2.*min_l.*(min_l+1).*(1+(-1).^(l(i)+l(ind_mn)))./(2.*m(i)+1)./(2.*n(i)+1);
    YY(i,ind_ln)=2.*min_m.*(min_m+1).*(1+(-1).^(m(i)+m(ind_ln)))./(2.*l(i)+1)./(2.*n(i)+1);
    ZZ(i,ind_lm)=2.*min_n.*(min_n+1).*(1+(-1).^(n(i)+n(ind_lm)))./(2.*l(i)+1)./(2.*m(i)+1);
    
    b(i,:)=(mu/4).*[(1/L(1))*(1-(-1)^l(i))*(m(i)==0)*(n(i)==0),(1/L(2))*(l(i)==0)*(1-(-1)^m(i))*(n(i)==0),(1/L(3))*(l(i)==0)*(m(i)==0)*(1-(-1)^n(i))];
end
A=(1/N).*Plmn'*Plmn;
B=(mu/32).*((1/L(1)^2).*XX+(1/L(2)^2).*YY+(1/L(3)^2).*ZZ);
a=(-(N_L^3/N).*(X_scaled-Y_scaled)*Plmn)';
XYZlmn=(A+B)\(a+b);

X_aligned=X_scaled+(1/N_L^3).*XYZlmn'*Plmn'; 
X_aligned=(X_aligned+1).*(((Max-Min)./2*ones(1,N)))+Min*ones(1,N);
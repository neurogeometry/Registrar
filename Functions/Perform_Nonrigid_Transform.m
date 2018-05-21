% This function performs non-rigid transformation (X -> X') based on Optimal_Nonrigid_Transform 
% X and X' are 3xN
% XYZlmn, N_L, Min, and Max are transformation parameters

function X_prime=Perform_Nonrigid_Transform(X,XYZlmn,N_L,Min,Max)

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
Plmn=zeros(N,N_L^3);
PP=P(X_scaled);
[l,m,n]=ind2sub([N_L,N_L,N_L],1:N_L^3);
for i=1:N_L^3
    Plmn(:,i)=PP(1,:,l(i)).* PP(2,:,m(i)).*PP(3,:,n(i));
end
X_prime=X_scaled+(1/N_L^3).*XYZlmn'*Plmn';
X_prime=(X_prime+1).*(((Max-Min)./2*ones(1,N)))+Min*ones(1,N);

% del_X2=sum((X_prime-X).^2,1);
% disp([mean(del_X2)^0.5, max(del_X2)^0.5])

% This function works with AM, AMlbl for branches, or AMlbl for trees. 
% r is fixed, while R is being optimized.
% Fitness = int(LoG(I)-alpha*(dR/ds)^2)ds

function R = Optimize_R(Orig,AMlbl,r,R,Rtypical,MaxIterations,alpha,betta,output)

if output==1
    disp('Optimization of R started.')
    format short g
    display(['        Iteration   ',   'I-cost   ',   '    R-cost   ',    ' Total Fitness'])
end

if betta>1
    betta=1;
end
if Rtypical<=2
    Rtypical=2;
end

epsilon=1; % added for stability
MinChange=betta*10^-3;
Rmin=1;
Rmax=5;

AM = spones(AMlbl+AMlbl');
N=length(AM);

[it,jt]=find(AM);
lll=AM; lll(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^0.5;
lllm1=AM; lllm1(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^(-0.5);
l=sum(lll,1)./2;
d2RcostdR2=2.*(diag(sum(lllm1,2))-lllm1);

if isempty(R) || sum(R)==0
    R=ones(N,1).*Rtypical;
end
R(R<Rmin)=Rmin;
R(R>Rmax)=Rmax;

W=ceil(2.5*Rmax)*2+1;
HW=(W-1)/2;
S=[W,W,W];
[xtemp,ytemp,ztemp]=ind2sub(S,1:prod(S));

%Orig = Multi_Scale_LoG(Orig,1,1,Rmax);
Orig=double(Orig);
Orig=Orig./max(Orig(:));
M=mean(Orig(:));
sizeIm=size(Orig);

Itemp=zeros(N,prod(S));
for i=1:N,
    xtemp1=xtemp+round(r(i,1))-HW-1;
    ytemp1=ytemp+round(r(i,2))-HW-1;
    ztemp1=ztemp+round(r(i,3))-HW-1;
    temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
    indIm=sub2ind_ASfast(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
    
    Im_S=ones(S).*M;
    Im_S(temp_ind)=Orig(indIm);
    Itemp(i,:)=Im_S(:)';
end

x=((r(:,1)-round(r(:,1)))*ones(1,prod(S))+HW+1-ones(N,1)*xtemp);
y=((r(:,2)-round(r(:,2)))*ones(1,prod(S))+HW+1-ones(N,1)*ytemp);
z=((r(:,3)-round(r(:,3)))*ones(1,prod(S))+HW+1-ones(N,1)*ztemp);
r2=(x.*x+y.*y+z.*z);

ChangeIcost=MinChange;
ChangeRcost=MinChange;
delR=zeros(N,1);
count=1;
while (ChangeIcost>=MinChange || ChangeRcost>=MinChange) && count<=MaxIterations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Intensity cost         
    r2R2=r2./repmat(R.^2,1,prod(S));
    r4R4=r2R2.*r2R2;
    r6R6=r2R2.*r4R4;
    
    LoG=exp(-r2R2).*(1-(2/3).*r2R2);
    d_LoG_dR=exp(-r2R2).*(-3+(16/3).*r2R2-(4/3).*r4R4);
    d2_LoG_dR2=exp(-r2R2).*(12-38.*r2R2+(64/3).*r4R4-(8/3).*r6R6);
    
    I_LoG=sum(LoG.*Itemp,2)./R.^3./pi^1.5;
    I_d_LoG_dR=sum(d_LoG_dR.*Itemp,2)./R.^4./pi^1.5;
    I_d2_LoG_dR2=sum(d2_LoG_dR2.*Itemp,2)./R.^5./pi^1.5; 
    
    Icost=l*I_LoG;
    dIcostdR=l'.*I_d_LoG_dR;
    d2IcostdR2=sparse(diag(l'.*I_d2_LoG_dR2-epsilon.*ones(N,1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % R cost
    RRR=AM; 
    RRR(AM>0)=(R(it)-R(jt));
    Rcostij=(lllm1.*RRR.^2);
    Rcost=sum(Rcostij(:))./2;
    dRcostdR=2.*sum(lllm1.*RRR,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fitness=Icost-alpha.*Rcost;
    if count==1 || Fitness>oldFitness
        if count>1
            ChangeIcost=abs((Icost-oldIcost)/Icost);
            ChangeRcost=abs((Rcost-oldRcost)/Rcost);
%             if ChangeIcost>0.1 || ChangeRcost>0.1
%                 error('The algorithm failed to converge. Decrease Optimization Step Size.')
%             end
        end
        oldIcost=Icost;
        oldRcost=Rcost;
        oldFitness=Fitness;
        Rold=R;
        delR=(d2IcostdR2-alpha.*d2RcostdR2)\(dIcostdR-alpha.*dRcostdR);
        R=R-betta.*delR;
        R(R<Rmin)=Rmin;
        R(R>Rmax)=Rmax;
    else
        ChangeIcost=abs((Icost-oldIcost)/Icost);
        ChangeRcost=abs((Rcost-oldRcost)/Rcost);
        Icost=oldIcost;
        Rcost=oldRcost;
        Fitness=oldFitness;
        betta=betta./1.5;
        R=Rold-betta.*delR;
        R(R<Rmin)=Rmin;
        R(R>Rmax)=Rmax;
    end
    if output==1
        disp(full([count, Icost, Rcost, Fitness]))
    end
    count=count+1;
end

if output==1
    disp('Optimization of R is complete.')
    if ChangeIcost>MinChange || ChangeRcost>MinChange
        disp('The algorithm did not converge to a solution.') 
        disp('You may increase Optimization Step Size and/or Maximum Number of Iterations and continue optimization.')
    end
end



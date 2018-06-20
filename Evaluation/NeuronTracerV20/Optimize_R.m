% This function works with AM, AMlbl for branches, or AMlbl for trees. 
% r is fixed, while R is being optimized.
% Fitness = int(LoG(I)-alpha*(dR/ds)^2)ds
% If Rtypical==0, R (if exists) will be used as starting configuration
% If Rtypical>0, R=Rtypical will be used as starting configuration

function R = Optimize_R(Orig,AMlbl,r,R,Rtypical,MaxIterations,alpha,betta,output,varargin)

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

epsilon=1; % added for stability
MinChange=betta*10^-3;
Rmin=0.5;
Rmax=30;

% if betta>1
%     betta=1;
% end

if Rtypical==0 % start with current R
    if isempty(R) || sum(R)==0
        Rtypical=2;
        R=ones(length(AMlbl),1).*Rtypical;
    end
else % start with Rtypical
    if Rtypical<=2
        Rtypical=2;
    elseif Rtypical>Rmax
        error(['This algorithm does not work for Typical Radius > ',num2str(Rmax),'. Consider reducing the image.'])
    end
    R=ones(length(AMlbl),1).*Rtypical;
end
R(R<Rmin)=Rmin;
R(R>Rmax)=Rmax;

AM = spones(AMlbl+AMlbl');

if output==1
    disp('Optimization of R started.')
    format short g
    display(['        Iteration   ',   'I-cost   ',   '    R-cost   ',    ' Total Fitness'])
end

[it,jt]=find(AM);
N=length(AM);
lll=AM; lll(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^0.5;
lllm1=AM; lllm1(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^(-0.5);
l=sum(lll,1)./2;
d2RcostdR2=2.*sparse(diag(sum(lllm1,2))-lllm1);
round_r=round(r);

Wmax=ceil(3*Rmax)*2+1;
NW=(Wmax+1)/2;
xtemp=cell(1,NW);
ytemp=cell(1,NW);
ztemp=cell(1,NW);
Wthr=41; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subsample=ones(1,NW);
for i=1:NW
    W=i*2-1;
    temp=1:W;
    if W>Wthr
        temp(rand(1,W)>(Wthr/W))=[];
        Subsample(i)=(length(temp)/W)^3;
    end
    [xtemp{i},ytemp{i},ztemp{i}]=ndgrid(temp,temp,temp);
end

Orig=double(Orig);
Orig=Orig./max(Orig(:));
M=mean(Orig(:));

sizeIm=size(Orig);
if length(sizeIm)==2
    sizeIm(3)=1;
end

ChangeIcost=MinChange;
ChangeRcost=MinChange;
delR=zeros(N,1);
I_LoG=zeros(N,1);
I_d_LoG_dR=zeros(N,1);
I_d2_LoG_dR2=zeros(N,1);
count=1;
while (ChangeIcost>=MinChange || ChangeRcost>=MinChange) && count<=MaxIterations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Intensity cost
    for i=1:N
        if myClass.getFlag()==1
            return;
        end
        
        W=ceil(3*R(i))*2+1;
        NW=(W+1)/2;
        
        xtemp1=xtemp{NW}(:)+round_r(i,1)-NW;
        ytemp1=ytemp{NW}(:)+round_r(i,2)-NW;
        ztemp1=ztemp{NW}(:)+round_r(i,3)-NW;
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=xtemp1(temp_ind)+(ytemp1(temp_ind)-1).*sizeIm(1)+(ztemp1(temp_ind)-1).*(sizeIm(1)*sizeIm(2));
        
        Itemp=ones(size(xtemp1)).*M;
        Itemp(temp_ind)=Orig(indIm);
        
        x=(r(i,1)-round_r(i,1)+NW-xtemp{NW}(:));
        y=(r(i,2)-round_r(i,2)+NW-ytemp{NW}(:));
        z=(r(i,3)-round_r(i,3)+NW-ztemp{NW}(:));
        r2R2=(x.*x+y.*y+z.*z).*(R(i)^-2);
        r4R4=r2R2.*r2R2;
        r6R6=r2R2.*r4R4;
        
        EXP=exp(-r2R2);
        LoG=EXP.*(1-(2/3).*r2R2);
        d_LoG_dR=EXP.*(-3+(16/3).*r2R2-(4/3).*r4R4);
        d2_LoG_dR2=EXP.*(12-38.*r2R2+(64/3).*r4R4-(8/3).*r6R6);
        
        I_LoG(i)=sum(LoG.*Itemp)/(R(i)^3*pi^1.5*Subsample(NW));
        I_d_LoG_dR(i)=sum(d_LoG_dR.*Itemp)/(R(i)^4*pi^1.5*Subsample(NW));
        I_d2_LoG_dR2(i)=sum(d2_LoG_dR2.*Itemp)/(R(i)^5*pi^1.5*Subsample(NW));
    end
       
    Icost=l*I_LoG;
    dIcostdR=l'.*I_d_LoG_dR;
    d2IcostdR2=diag(sparse(l'.*I_d2_LoG_dR2-epsilon.*ones(N,1)));
    
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
            ChangeIcost=abs((Icost-oldIcost)/oldIcost);
            ChangeRcost=abs((Rcost-oldRcost)/oldRcost);
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
        ChangeIcost=abs((Icost-oldIcost)/oldIcost);
        ChangeRcost=abs((Rcost-oldRcost)/oldRcost);
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

R=Rold;

if output==1
    disp('Optimization of R is complete.')
    if ChangeIcost>MinChange || ChangeRcost>MinChange
        disp('The algorithm did not converge to solution with default precision.') 
        disp('Consider increasing Optimization Step Size and/or Maximum Number of Iterations.')
    end
end



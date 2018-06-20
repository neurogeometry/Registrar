% This function works with AM, AMlbl for branches, or AMlbl for trees. 
% Trees are be optimized separately.
% Branch positions (r) and calibers (R) are optimized simultaneously  
% Branch and end points can be fixed or optimized:
% Optimize_bps = 1,0 optimize branch points.
% Optimize_tps = 1,0 optimize terminal (start, end) points.
% AMlbl in the output is labled for trees 
% This version of the code is normalized for ppm as in the paper
% If Rtypical==0, R (if exists) will be used as starting configuration
% If Rtypical>0, R=Rtypical will be used as starting configuration

function [AMlbl r R I_F]=Optimize_Trace_R(Orig,AMlbl,r,R,Rtypical,Optimize_bps,Optimize_tps,pointsperum,MaxIterations,alpha_r,alpha_R,betta,adjustPPM,output,varargin)

if length(varargin)==1
    myClass=varargin{1};
else
    myClass.getFlag=NaN;
end

epsilon=1; % added for stability
MinChange=betta*10^-6;
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

% adjust vertex density
if adjustPPM==1
    [AM r R] = AdjustPPM(AMlbl,r,R,pointsperum);
    AM=spones(AM);
else
    AM=spones(AMlbl);
end

if output==1
    disp('F1 trace optimization started.')
    format short g
    display(['        Iteration   ',   'I-cost   ',   '    L-cost   ',   '    R-cost   ',    ' Total Fitness'])
end

[it,jt]=find(AM);
N=length(AM);
B=sparse(diag(sum(AM,2))-AM);
B4=2.*pointsperum.*blkdiag(alpha_r.*B,alpha_r.*B,alpha_r.*B,alpha_R.*B);

Wmax=ceil(3*Rmax)*2+1;
NW=(Wmax+1)/2;
xtemp=cell(1,NW);
ytemp=xtemp;
ztemp=xtemp;
Wthr=41; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subsample=ones(1,NW);
for i=1:NW
    W=i*2-1;
    temp=-(i-1):(i-1);
    if W>Wthr
        temp(rand(1,W)>(Wthr/W))=[];
        Subsample(i)=(length(temp)/W)^3;
    end
    [xtemp{i},ytemp{i},ztemp{i}]=ndgrid(temp,temp,temp);
%     xtemp{i}=repmat(temp',[1,length(temp),length(temp)]);
%     ytemp{i}=permute(xtemp{i},[2,1,3]);
%     ztemp{i}=permute(xtemp{i},[3,2,1]);
end

sizeIm=size(Orig);
if length(sizeIm)==2
    sizeIm(3)=1;
end

Orig=double(Orig(:));
Orig=Orig./max(Orig);
M=mean(Orig);

tps=(sum(AM)==1);
bps=(sum(AM)>2);
move=true(N,1);
if Optimize_tps==0
    move(tps)=false;
end
if Optimize_bps==0
    move(bps)=false;
end

ChangeIcost=MinChange;
ChangeLcost=MinChange;
ChangeRcost=MinChange;
del_r=zeros(size(r));
del_R=zeros(N,1);
I_Q=zeros(N,1); I_Qx=zeros(N,1); I_Qy=zeros(N,1); I_Qz=zeros(N,1);
I_Pxx=zeros(N,1); I_Pxy=zeros(N,1); I_Pxz=zeros(N,1); I_Pyy=zeros(N,1); I_Pyz=zeros(N,1); I_Pzz=zeros(N,1);
I_F=zeros(N,1); I_J=zeros(N,1); I_H=zeros(N,1);
I_Kx=zeros(N,1); I_Ky=zeros(N,1); I_Kz=zeros(N,1);
count=1;
while (ChangeIcost>=MinChange || ChangeLcost>=MinChange || ChangeRcost>=MinChange) && count<=MaxIterations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Intensity cost, LoG filter
    for i=1:N
        if myClass.getFlag()==1
            return;
        end
    
        W=ceil(3*R(i))*2+1;
        NW=(W+1)/2;
        
        xtemp1=xtemp{NW}(:)+round(r(i,1));
        ytemp1=ytemp{NW}(:)+round(r(i,2));
        ztemp1=ztemp{NW}(:)+round(r(i,3));
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=xtemp1(temp_ind)+(ytemp1(temp_ind)-1).*sizeIm(1)+(ztemp1(temp_ind)-1).*(sizeIm(1)*sizeIm(2));
        
        Itemp=ones(size(xtemp1)).*M;
        Itemp(temp_ind)=Orig(indIm);
        
        temp1=1/R(i);
        xR=(r(i,1)-xtemp1)*temp1;
        yR=(r(i,2)-ytemp1)*temp1;
        zR=(r(i,3)-ztemp1)*temp1;
        r2R2=(xR.*xR+yR.*yR+zR.*zR);
        r4R4=r2R2.*r2R2;
        r6R6=r2R2.*r4R4;
        
        EXP=exp(-r2R2).*Itemp;
        F=EXP.*(1-(2/3).*r2R2);
        Q=EXP.*(-10/3+(4/3).*r2R2);
        P=EXP.*(28/3-(8/3).*r2R2);
        H=EXP.*(-3+(16/3).*r2R2-(4/3).*r4R4);
        J=EXP.*(12-38.*r2R2+(64/3).*r4R4-(8/3).*r6R6);
        K=EXP.*(50/3-16.*r2R2+(8/3).*r4R4);
        
        temp3=1/(R(i)^3*pi^1.5*pointsperum*Subsample(NW));
        temp4=1/(R(i)^4*pi^1.5*pointsperum*Subsample(NW));
        temp5=1/(R(i)^5*pi^1.5*pointsperum*Subsample(NW));
        
        I_Q(i)=sum(Q)*temp5;
        I_Qx(i)=sum(Q.*xR)*temp4;
        I_Qy(i)=sum(Q.*yR)*temp4;
        I_Qz(i)=sum(Q.*zR)*temp4;
        
        Px=P.*xR;
        Py=P.*yR;
        Pz=P.*zR;
        I_Pxx(i)=sum(Px.*xR)*temp5;
        I_Pxy(i)=sum(Px.*yR)*temp5;
        I_Pxz(i)=sum(Px.*zR)*temp5;
        I_Pyy(i)=sum(Py.*yR)*temp5;
        I_Pyz(i)=sum(Py.*zR)*temp5;
        I_Pzz(i)=sum(Pz.*zR)*temp5;
        
        I_F(i)=sum(F)*temp3;
        I_H(i)=sum(H)*temp4;
        I_J(i)=sum(J)*temp5;
        
        I_Kx(i)=sum(K.*xR)*temp5;
        I_Ky(i)=sum(K.*yR)*temp5;
        I_Kz(i)=sum(K.*zR)*temp5;
    end
    
    Icost=sum(I_F);
    dFcostdrR=[I_Qx;I_Qy;I_Qz;I_H]-B4*[r(:);R];
    d2IcostdRr2=[diag(sparse(I_Pxx+I_Q)),diag(sparse(I_Pxy)),diag(sparse(I_Pxz)),diag(sparse(I_Kx));
            diag(sparse(I_Pxy)),diag(sparse(I_Pyy+I_Q)),diag(sparse(I_Pyz)),diag(sparse(I_Ky));
            diag(sparse(I_Pxz)),diag(sparse(I_Pyz)),diag(sparse(I_Pzz+I_Q)),diag(sparse(I_Kz));
            diag(sparse(I_Kx)),diag(sparse(I_Ky)),diag(sparse(I_Kz)),diag(sparse(I_J))]-B4-epsilon.*diag(sparse(ones(4*N,1)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % L cost
    Lcost=sum(sum((r(it,:)-r(jt,:)).^2,2))/2*pointsperum;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % R cost
    Rcost=sum((R(it)-R(jt)).^2)/2*pointsperum;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fitness=Icost-alpha_r.*Lcost-alpha_R.*Rcost;
    if count==1 || Fitness>oldFitness
        if count>1
            ChangeIcost=abs((Icost-oldIcost)/oldIcost);
            ChangeLcost=abs((Lcost-oldLcost)/oldLcost);
            ChangeRcost=abs((Rcost-oldRcost)/oldRcost);
            if ChangeIcost>0.5 || ChangeLcost>0.5
                warning('Trace may be unstable. Decrease Optimization Step Size and/or increase Trace Stiffness parameters.')
            end
        end
        oldIcost=Icost;
        oldLcost=Lcost;
        oldRcost=Rcost;
        oldFitness=Fitness;
        r_old=r;
        R_old=R;
        
        del_rR=d2IcostdRr2\dFcostdrR;
        del_r(:)=del_rR(1:3*N);
        del_R=del_rR(3*N+1:end);
        
        r(move,:)=r(move,:)-betta.*del_r(move,:);
        r(r(:,1)<0.5,1)=0.5; r(r(:,1)>sizeIm(1)+0.5,1)=sizeIm(1)+0.5; % r ranges from 0.5 to sizeIm+0.5 in Matlab and [0 sizeIm] in Java and SWC
        r(r(:,2)<0.5,2)=0.5; r(r(:,2)>sizeIm(2)+0.5,2)=sizeIm(2)+0.5;
        r(r(:,3)<0.5,3)=0.5; r(r(:,3)>sizeIm(3)+0.5,3)=sizeIm(3)+0.5;
                
        R=R-betta.*del_R;
        R(R<Rmin)=Rmin;
        R(R>Rmax)=Rmax;
    else
        ChangeIcost=abs((Icost-oldIcost)/oldIcost);
        ChangeLcost=abs((Lcost-oldLcost)/oldLcost);
        ChangeRcost=abs((Rcost-oldRcost)/oldRcost);
        Icost=oldIcost;
        Lcost=oldLcost;
        Rcost=oldRcost;
        Fitness=oldFitness;
        betta=betta./1.5;
        r(move,:)=r_old(move,:)-betta.*del_r(move,:);
        r(r(:,1)<0.5,1)=0.5; r(r(:,1)>sizeIm(1)+0.5,1)=sizeIm(1)+0.5;
        r(r(:,2)<0.5,2)=0.5; r(r(:,2)>sizeIm(2)+0.5,2)=sizeIm(2)+0.5;
        r(r(:,3)<0.5,3)=0.5; r(r(:,3)>sizeIm(3)+0.5,3)=sizeIm(3)+0.5;
        R=R_old-betta.*del_R;
        R(R<Rmin)=Rmin;
        R(R>Rmax)=Rmax;
    end
 
    if output==1
        disp(full([count, Icost, Lcost, Rcost, Fitness]))
    end
    count=count+1;
end

r=r_old;
R=R_old;
AMlbl = LabelTreesAM(AM);

if output==1
    disp('F1 trace optimization is complete.')
    if ChangeIcost>MinChange || ChangeLcost>MinChange || ChangeRcost>MinChange
        disp('The algorithm did not converge to solution with default precision.') 
        disp('Consider increasing Optimization Step Size and/or Maximum Number of Iterations.')
    end   
end

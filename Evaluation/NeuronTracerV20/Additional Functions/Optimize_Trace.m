% This function works with AM, AMlbl for branches, or AMlbl for trees. 
% Trees are be optimized separately.
% Branch and end points can be fixed or optimized:
% Optimize_bps = 1,0 optimize branch points.
% Optimize_tps = 1,0 optimize terminal (start, end) points.
% Rtypical is the radius of branches which remains fixed.
% AMlbl in the output is labled for trees 
% This version of the code is normalized for ppm as in the paper

function [AMlbl r I_snake]=Optimize_Trace(Orig,AMlbl,r,Rtypical,Optimize_bps,Optimize_tps,pointsperum,MaxIterations,alpha,betta,adjustPPM,output)

epsilon=1; % added for stability
MinChange=betta*10^-3;
Rmin=1;
Rmax=30;

if betta>1
    betta=1;
end

if Rtypical<=Rmin
    Rtypical=Rmin;
elseif Rtypical>Rmax
    error(['This algorithm does not work for Typical Radius > ',num2str(Rmax),'. Consider reducing the image.'])
end
R=ones(length(AMlbl),1).*Rtypical;

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
    display(['        Iteration   ',   'I-cost   ',   '    L-cost   ',    ' Total Fitness'])
end

[it,jt]=find(AM);
N=length(AM);
d2Lcostdr2=pointsperum.*2.*(diag(sum(AM,2))-AM);

W=ceil(3*Rtypical)*2+1;
HW=(W-1)/2;
S=[W,W,W];
[xtemp,ytemp,ztemp]=ind2sub(S,1:prod(S));

Orig=double(Orig);
Orig=Orig./max(Orig(:));
M=mean(Orig(:));
sizeIm=size(Orig);

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
del_r=zeros(size(r));
count=1;
while (ChangeIcost>=MinChange || ChangeLcost>=MinChange) && count<=MaxIterations
    Itemp=zeros(N,prod(S)); 
    for i=1:N
        xtemp1=xtemp+round(r(i,1))-HW-1;
        ytemp1=ytemp+round(r(i,2))-HW-1;
        ztemp1=ztemp+round(r(i,3))-HW-1;
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=xtemp1(temp_ind)+(ytemp1(temp_ind)-1).*sizeIm(1)+(ztemp1(temp_ind)-1).*(sizeIm(1)*sizeIm(2));
        
        Im_S=ones(S).*M; 
        Im_S(temp_ind)=Orig(indIm);
        Itemp(i,:)=Im_S(:)'; 
    end
    
    xR=(repmat(r(:,1)-round(r(:,1))+HW+1,1,prod(S))-repmat(xtemp,N,1))./repmat(R,1,prod(S));
    yR=(repmat(r(:,2)-round(r(:,2))+HW+1,1,prod(S))-repmat(ytemp,N,1))./repmat(R,1,prod(S));
    zR=(repmat(r(:,3)-round(r(:,3))+HW+1,1,prod(S))-repmat(ztemp,N,1))./repmat(R,1,prod(S));
    r2R2=(xR.*xR+yR.*yR+zR.*zR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Intensity cost
    
    %%%%%%% Gaussian filter
%     F=exp(-r2R2).*Itemp;
%     Q=EXP.*(-10/3+(4/3).*r2R2);
%     P=EXP.*(28/3-(8/3).*r2R2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% LoG filter
    EXP=exp(-r2R2).*Itemp;
    F=EXP.*(1-(2/3).*r2R2);
    Q=EXP.*(-10/3+(4/3).*r2R2);
    P=EXP.*(28/3-(8/3).*r2R2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    I_Q=sum(Q,2)./R.^5./pi^1.5/pointsperum-epsilon*sparse(ones(N,1));
    I_Qx=sum(Q.*xR,2)./R.^4./pi^1.5/pointsperum;
    I_Qy=sum(Q.*yR,2)./R.^4./pi^1.5/pointsperum;
    I_Qz=sum(Q.*zR,2)./R.^4./pi^1.5/pointsperum;
    
    Px=P.*xR;
    Py=P.*yR;
    Pz=P.*zR;
    I_Pxx=sum(Px.*xR,2)./R.^5./pi^1.5/pointsperum;
    I_Pxy=sum(Px.*yR,2)./R.^5./pi^1.5/pointsperum;
    I_Pxz=sum(Px.*zR,2)./R.^5./pi^1.5/pointsperum;
    I_Pyy=sum(Py.*yR,2)./R.^5./pi^1.5/pointsperum;
    I_Pyz=sum(Py.*zR,2)./R.^5./pi^1.5/pointsperum;
    I_Pzz=sum(Pz.*zR,2)./R.^5./pi^1.5/pointsperum;
    
    I_F=sum(F,2)./R.^3./pi^1.5;
    Icost=sum(I_F)/pointsperum;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % L cost
    lll2=AM; lll2(AM>0)=sum((r(it,:)-r(jt,:)).^2,2);
    Lcost=sum(lll2(:))/2*pointsperum;
    dLcostdr=(2*pointsperum).*(diag(sum(AM,2))-AM)*r;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fitness=Icost-alpha.*Lcost;
    if count==1 || Fitness>oldFitness
        if count>1
            ChangeIcost=abs((Icost-oldIcost)/Icost);
            ChangeLcost=abs((Lcost-oldLcost)/Lcost);
            if ChangeIcost>0.5 || ChangeLcost>0.5
                error('Trace is unstable. Decrease Optimization Step Size and/or Trace Stiffness.')
            end
        end
        oldIcost=Icost;
        oldLcost=Lcost;
        oldFitness=Fitness;
        r_old=r;
        del_r(:)=([diag(sparse(I_Pxx+I_Q))-alpha.*d2Lcostdr2,diag(sparse(I_Pxy)),diag(sparse(I_Pxz));
            diag(sparse(I_Pxy)),diag(sparse(I_Pyy+I_Q))-alpha.*d2Lcostdr2,diag(sparse(I_Pyz));
            diag(sparse(I_Pxz)),diag(sparse(I_Pyz)),diag(sparse(I_Pzz+I_Q))-alpha.*d2Lcostdr2])\([I_Qx;I_Qy;I_Qz]-alpha.*dLcostdr(:));
        r(move,:)=r(move,:)-betta.*del_r(move,:);
        r(r(:,1)<0.5,1)=0.5; r(r(:,1)>sizeIm(1)+0.5,1)=sizeIm(1)+0.5;
        r(r(:,2)<0.5,2)=0.5; r(r(:,2)>sizeIm(2)+0.5,2)=sizeIm(2)+0.5;
        r(r(:,3)<0.5,3)=0.5; r(r(:,3)>sizeIm(3)+0.5,3)=sizeIm(3)+0.5;
    else
        ChangeIcost=abs((Icost-oldIcost)/Icost);
        ChangeLcost=abs((Lcost-oldLcost)/Lcost);
        Icost=oldIcost;
        Lcost=oldLcost;
        Fitness=oldFitness;
        betta=betta./1.5;
        r(move,:)=r_old(move,:)-betta.*del_r(move,:);
        r(r(:,1)<0.5,1)=0.5; r(r(:,1)>sizeIm(1)+0.5,1)=sizeIm(1)+0.5;
        r(r(:,2)<0.5,2)=0.5; r(r(:,2)>sizeIm(2)+0.5,2)=sizeIm(2)+0.5;
        r(r(:,3)<0.5,3)=0.5; r(r(:,3)>sizeIm(3)+0.5,3)=sizeIm(3)+0.5;
    end
 
    if output==1
        disp(full([count, Icost, Lcost, Fitness]))
    end
    count=count+1;
end

r=r_old;
AMlbl = LabelTreesAM(AM);
I_snake=I_F./pointsperum;

if output==1
    disp('F1 trace optimization is complete.')
    if ChangeIcost>MinChange || ChangeLcost>MinChange
        disp('The algorithm did not converge to solution with default precision.') 
        disp('Consider increasing Optimization Step Size and/or Maximum Number of Iterations.')
    end   
end

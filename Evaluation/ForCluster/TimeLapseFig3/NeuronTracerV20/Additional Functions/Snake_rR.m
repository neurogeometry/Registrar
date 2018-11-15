% This version works with AM or AMlbl. Trees are optimized separately.
% r and R are optimized at the same time
% Branch and end points can also be optimized
% Optimize_bps = 1,0 optimize branch points.
% Optimize_tps = 1,0 optimize terminal (start, end) points.
% Multiple_trees = 1,0 optimize labeled trees separately.
% the cost is normalized for ppm

function [AMlbl r R I_snake]=Snake_rR(Orig,AMlbl,r,R,R_min,R_max,Optimize_bps,Optimize_tps,Multiple_trees,pointsperum,Nstep,alpha,betta,gamma,output)

if output==1
    disp('Trace optimization started.')
    format short g
    display(['      Iteration   ',   '# vertices   ',   'Trace length   ',    '<Intensity>'])
end

R_step=0.1;
if isempty(R)
    R=Find_R_fast(Orig,r,R_min,R_step,R_max);
end

AMlbl=max(AMlbl,AMlbl');
rem_ind=(sum(AMlbl,1)==0);
AMlbl(:,rem_ind)=[];
AMlbl(rem_ind,:)=[];
r(rem_ind,:)=[];
R(rem_ind)=[];
AM = spones(AMlbl);

if Multiple_trees==1
    % disconnect different label trees (if connected)
    bp=find(sum(AM)>2);
    for i=1:length(bp)
        bp_labels=nonzeros(AMlbl(:,bp(i)));
        ubp_labels=unique(bp_labels);
        if length(ubp_labels)>1
            r=[r;ones(length(ubp_labels)-1,1)*r(bp(i),:)];
            R=[R;R(bp(i)).*ones(length(ubp_labels)-1,1)];
            a=length(AM);
            AM(a+length(ubp_labels)-1,a+length(ubp_labels)-1)=0;
            for j=2:length(ubp_labels)
                temp=(AMlbl(:,bp(i))==ubp_labels(j));
                AM(temp,bp(i))=0;
                AM(bp(i),temp)=0;
                AM(temp,a+j-1)=1;
                AM(a+j-1,temp)=1;
            end
        end
    end
end

[it,jt]=find(AM);
lll=AM; lll(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^0.5;
L0=sum(lll(:))/2;

[AMlbl r R] = AdjustPPM(AM,r,R,pointsperum);
AM = spones(AMlbl);

Orig = Multi_Scale_Gaussian(Orig,mean(R),1,1);
Orig=double(Orig);
Orig=Orig./max(Orig(:));
sizeIm=size(Orig);

Imx=zeros(sizeIm); Imx(2:end-1,:,:)=(Orig(3:end,:,:)-Orig(1:end-2,:,:))./2;
Imy=zeros(sizeIm); Imy(:,2:end-1,:)=(Orig(:,3:end,:)-Orig(:,1:end-2,:))./2;
Imz=zeros(sizeIm); Imz(:,:,2:end-1)=(Orig(:,:,3:end)-Orig(:,:,1:end-2))./2;

for f1=1:Nstep
    % resegment
%     [AM rt Rt] = Divide_Segments(AM,rt,Rt',pointsperum);
%     [AM rt Rt] = Merge_Segments(AM,rt,Rt',pointsperum);
    N=length(AM);
    [it,jt]=find(AM);
    lll=AM; lll(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^0.5;
    lllm1=AM; lllm1(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^(-0.5);
    L=sum(lll(:))/2;
    if L/L0>3
        error('Trace is unstable. Reduce Trace Stiffness (alpha) or Optimization Step Size (beta).')
    end
    
    RRR=AM; RRR(AM>0)=(R(it)-R(jt));
    temp=sum(AM);
    intermp=find(temp==2);
    endp=find(temp==1);
    bp=find(temp>2);
         
    if Optimize_bps==1
        
    end
    
    if Optimize_tps==1
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Intensity gradients
    I_snake=interp3(Orig,r(:,2),r(:,1),r(:,3),'linear',0);
    dIcostdr=zeros(N,3);
    dIcostdr(:,1)=interp3(Imx,r(:,2),r(:,1),r(:,3),'linear',0);
    dIcostdr(:,2)=interp3(Imy,r(:,2),r(:,1),r(:,3),'linear',0);
    dIcostdr(:,3)=interp3(Imz,r(:,2),r(:,1),r(:,3),'linear',0);
    dIcostdR=zeros(N,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Lennard-Jones potential gradients
    lmin=1/pointsperum;
    Ucostij=(lllm1.*lmin).^(12)-2.*(lllm1.*lmin).^(6);
    Ucost=sum(Ucostij(:))./2;
    temp=-(12./lmin^2).*((lllm1.*lmin).^(14)-(lllm1.*lmin).^(8));
    dUcostdr=(diag(sum(temp,1))-temp)*r;
    dUcostdR=zeros(N,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rcost gradients
    Rcostij=(lllm1.*RRR.^2);
    Rcost=sum(Rcostij(:))./2;
    temp=(lllm1.^3.*RRR.^2);
    dRcostdr=-(diag(sum(temp,1))-temp)*r;
    dRcostdR=full(2.*sum(lllm1.*RRR))';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if output==1
        DgIav=sum(diag(sum(lll))/sum(lll(:))*I_snake);
        disp(full([f1, N, L, DgIav]))
    end
    
    r=r+0.01.*(dIcostdr./pointsperum-0.01*dUcostdr-0.01*dRcostdr);
    R=R+0.01.*(dIcostdR./pointsperum-1*dUcostdR+1*dRcostdR);
    
    if Optimize_bps==1
        
    end
    if Optimize_tps==1
        
    end
    
    r(r(:,1)<0,1)=0; r(r(:,2)<0,2)=0; r(r(:,3)<0,3)=0;
    r(r(:,1)>sizeIm(1),1)=sizeIm(1); r(r(:,2)>sizeIm(2),2)=sizeIm(2); r(r(:,3)>sizeIm(3),3)=sizeIm(3);
    R(R<0)=R_min;

end

if Multiple_trees==0
    AMlbl = LabelBranchesAM(AM>0);
elseif Multiple_trees==1
    AMlbl = LabelTreesAM(AM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds the gradient of intensity at locations r
function [gradI,meanI]=grad_Im(Orig,rt,sigma)

sizeIm=size(Orig);
Orig=double(Orig);
Orig=Orig./max(Orig(:));

W=ceil(3*sigma)*2+1;
HW=(W-1)/2;
prodW=prod(W);
[xtemp,ytemp,ztemp]=ind2sub(W,1:prodW);
N=size(rt,2);

gradI=zeros(3,N); 
meanI=zeros(1,N);

Fxtemp=zeros(prodW,N); Fytemp=zeros(prodW,N); Fztemp=zeros(prodW,N);
Itemp=zeros(prodW,N);

Gx=exp(-((rt(1,:)-round(rt(1,:)))'*ones(1,prodS)+HW+1-ones(N,1)*xtemp).^2./2./sigma^2);
Gy=exp(-((rt(2,:)-round(rt(2,:)))'*ones(1,prodS)+HW+1-ones(N,1)*ytemp).^2./2./sigma^2);
Gz=exp(-((rt(3,:)-round(rt(3,:)))'*ones(1,prodS)+HW+1-ones(N,1)*ztemp).^2./2./sigma^2);
G=Gx.*Gy.*Gz;
G=G./(sum(G,2)*ones(1,prodS));

for i=1:N,
    xtemp1=xtemp+round(rt(1,i))-HW-1;
    ytemp1=ytemp+round(rt(2,i))-HW-1;
    ztemp1=ztemp+round(rt(3,i))-HW-1;
    temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
    indIm=sub2ind_ASfast(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
    
    Im_S=zeros(S);
    Im_S(temp_ind)=Orig(indIm);
    
    Fxx=Im_S;
    Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2;
    Fxtemp(:,i)=Fxx(:);
    
    Fyy=Im_S;
    Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2;
    Fytemp(:,i)=Fyy(:);
    
    Fzz=Im_S;
    Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2;
    Fztemp(:,i)=Fzz(:);
    
    Itemp(:,i)=Im_S(:);
end

gradI(1,intermp)=sum(G'.*Fxtemp);
gradI(2,intermp)=sum(G'.*Fytemp);
gradI(3,intermp)=sum(G'.*Fztemp);
I_snake(intermp)=sum(G'.*Itemp);

if output==1
    disp('Trace optimization is complete.')
end


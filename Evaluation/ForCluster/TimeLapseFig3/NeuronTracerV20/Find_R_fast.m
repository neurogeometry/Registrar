% This function finds best R for all verticies in the trace
% R is defined as sqrt(2) times the SD of the best LoG filter
% This version rounds vertex positions

function [R,I] = Find_R_fast(Orig,r,R_min,R_step,R_max)

%disp('Find R fast process started.')

if R_min<0.1
    R_min=0.1;
end

s=(R_min:R_step:R_max)';

if isempty(s)
    R=zeros(size(r(:,1)));
    I=zeros(size(r(:,1)))';
else
    r=round(r); % !!!!!!!!!
    del=0.01*R_min;
    N=size(r,1);
    
    sizeIm=size(Orig);
    if length(sizeIm)==2
        sizeIm=[sizeIm,1];
    end
    
    W=ceil(2*max(s))*2+1;
    HW=(W-1)/2;
    S=[W,W,W];
    [xtemp,ytemp,ztemp]=ind2sub(S,1:prod(S));
    
    Itemp=zeros(prod(S),N);
    for i=1:N,
        xtemp1=xtemp+r(i,1)-HW-1;
        ytemp1=ytemp+r(i,2)-HW-1;
        ztemp1=ztemp+r(i,3)-HW-1;
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=sub2ind_ASfast(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
        
        Im_S=zeros(S);
        Im_S(temp_ind)=double(Orig(indIm));
        
        Itemp(:,i)=Im_S(:);
    end
    
    r2=(HW+1-xtemp).^2+(HW+1-ytemp).^2+(HW+1-ztemp).^2;
    Gm=exp(-(1./(s-del).^2)*r2).*((1./(s-del).^3)*ones(1,prod(S)));
    Gm=Gm./(sum(Gm,2)*ones(1,prod(S)));
    Gp=exp(-(1./(s+del).^2)*r2).*((1./(s+del).^3)*ones(1,prod(S)));
    Gp=Gp./(sum(Gp,2)*ones(1,prod(S)));
    LoG=(Gm-Gp)./(2*del).*(s*ones(1,prod(S))); % correct up to a numerical factor
    
    I_LoG=LoG*Itemp;
    if length(s)>1
        [I,ind]=max(I_LoG);
        R=s(ind);
    else
        I=I_LoG;
        R=s.*ones(size(r(:,1)));
    end
end

%disp('Find R fast process is complete.')


% This function finds best R for all verticies in the trace
% R is defined as sqrt(2) times the SD of the best LoG filter

function R = Find_R(Orig,r,R_min,R_step,R_max)

disp('Find R process started.')

N=size(r,1);
sizeIm=size(Orig);

if R_min<0.1
    R_min=0.1;
end

s=(R_min:R_step:R_max);

if isempty(s)
    R=zeros(N,1);
else
    del=0.01*R_min;
    
    W=ceil(2*max(s))*2+1;
    HW=(W-1)/2;
    S=[W,W,W];
    [xtemp,ytemp,ztemp]=ind2sub(S,1:prod(S));
    
    Itemp=zeros(prod(S),N);
    for i=1:N,
        xtemp1=xtemp+round(r(i,1))-HW-1;
        ytemp1=ytemp+round(r(i,2))-HW-1;
        ztemp1=ztemp+round(r(i,3))-HW-1;
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=sub2ind_ASfast(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
        
        Im_S=zeros(S);
        Im_S(temp_ind)=double(Orig(indIm));
        
        Itemp(:,i)=Im_S(:);
    end
    
    x2=((r(:,1)-round(r(:,1)))*ones(1,prod(S))+HW+1-ones(N,1)*xtemp).^2;
    y2=((r(:,2)-round(r(:,2)))*ones(1,prod(S))+HW+1-ones(N,1)*ytemp).^2;
    z2=((r(:,3)-round(r(:,3)))*ones(1,prod(S))+HW+1-ones(N,1)*ztemp).^2;
    r2=(x2+y2+z2)';
    
    R=zeros(N,1);
    I_LoG_max=zeros(1,N);
    for i=1:length(s)
        %LoG=exp(-r2./s(i)^2).*(1-r2./(1.5*s(i)^2))./s(i)^3;
        
        Gm=exp(-r2./(s(i)-del)^2)./(s(i)-del)^3;
        Gm=Gm./(ones(prod(S),1)*sum(Gm,1));
        Gp=exp(-r2./(s(i)+del)^2)./(s(i)+del)^3;
        Gp=Gp./(ones(prod(S),1)*sum(Gp,1));

        LoG=(Gm-Gp)./(2*del).*s(i); % correct up to a numerical factor
        
        I_LoG_temp=sum(LoG.*Itemp);
        ind=(I_LoG_temp>I_LoG_max);
        I_LoG_max(ind)=I_LoG_temp(ind);
        R(ind)=s(i);
    end
end

disp('Find R process is complete.')


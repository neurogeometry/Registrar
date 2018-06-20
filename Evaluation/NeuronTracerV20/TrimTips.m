% This function finds terminal branches and trims them based on
% variance in intensity and contrast

function [AMlbl,r]=TrimTips(AMlbl,r,Orig,stepbackdist,R_min,R_step,R_max)

pad=10;
Orig=double(Orig);
sizeIm=size(Orig);
if length(sizeIm)==2 
    sizeIm=[sizeIm,1];
end

r_ind=sub2ind(sizeIm,r(:,1),r(:,2),r(:,3));
tips=StepBackTips(AMlbl,stepbackdist);
[~, I] = Find_R_fast(Orig,r,R_min,R_step,R_max);

%background=mean(Orig(:));
%foreground=mean(I);

deletepoints=false(size(I));
for i=1:length(tips)
    n=length(tips{i});
    Imval=I(tips{i});
    foreground=mean(Imval);
    %     cumsumIm=cumsum(Imval);
    %     cumsumImsq=cumsum(Imval.^2);
    %
    %     %calculating mean of squares
    %     LcumsumImsq=cumsumImsq(1:end-1);
    %     RcumsumImsq=cumsumImsq(end)-cumsumImsq(1:end-1);
    %
    %     %calculating square of means
    %     LcumsumIm=cumsumIm(1:end-1)./(1:n-1);
    %     RcumsumIm=(cumsumIm(end)-cumsumIm(1:end-1))./(n-1:-1:1);
    %
    %     %Normalized variances:
    %     Lvar=LcumsumImsq./(1:n-1)-LcumsumIm.^2;
    %     Rvar=RcumsumImsq./(n-1:-1:1)-RcumsumIm.^2;
    %     variance=[var(Imval,1),(Lvar.*(1:n-1)+Rvar.*(n-1:-1:1))/n];
    %
    %     contrast=[mean(Imval),RcumsumIm-LcumsumIm];
    %     contrast(contrast<0)=0;
    %
    %     cost=variance.^.5./contrast; %./(length(contrast):-1:1)'.^.5;
    
    [GDx,GDy,GDz]=ind2sub(sizeIm,tips{i});
    min_x=min(GDx); max_x=max(GDx);
    min_y=min(GDy); max_y=max(GDy);
    min_z=min(GDz); max_z=max(GDz);
    
    Xmin=max(min_x-1-pad,0)+1;
    Xmax=min(max_x-sizeIm(1)+pad,0)+sizeIm(1);
    Ymin=max(min_y-1-pad,0)+1;
    Ymax=min(max_y-sizeIm(2)+pad,0)+sizeIm(2);
    Zmin=max(min_z-1-pad,0)+1;
    Zmax=min(max_z-sizeIm(3)+pad,0)+sizeIm(3);
    
    Im1=Orig(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
    background=mean(Im1(:));
    
    Fsum=cumsum((foreground-Imval(n:-1:1)).^2);
    Bsum=cumsum((Imval-background).^2);
    cost=Fsum(n:-1:1)+[0,Bsum(1:end-1)];
    
%         figure(10),plot(Imval,'.-'),hold on
%         plot([1,length(Imval)],[background,background],'-')
%         plot([1,length(Imval)],[foreground,foreground],'-')
%         plot(cost./max(cost(~isinf(cost))).*foreground,'r-')
%         plot(cost./max(cost).*foreground,'g-')
%         hold off
    
    [~,cutoffind]=min(cost);
    
    if cutoffind>1
        deletepoints(tips{i}(1:cutoffind-1))=true;
    end
end

AMlbl(deletepoints,:)=[];
AMlbl(:,deletepoints)=[];
r_ind(deletepoints)=[];
[rx,ry,rz]=ind2sub(sizeIm,r_ind);
r=[rx,ry,rz];

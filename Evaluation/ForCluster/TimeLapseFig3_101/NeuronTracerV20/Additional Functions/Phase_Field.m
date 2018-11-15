% This function preprocesses the stack according to the phase-field method

function Im = Phase_Field(Im,T,alpha,h)

cl=class(Im);
Im=double(Im); 
maxIm=max(Im(:));
Im=Im./maxIm;
sizeIm=size(Im);
step_size=0.1;
ind_fixed=(Im(:)>0.8*maxIm);

Fgr=zeros(sizeIm);
Ff=zeros(sizeIm);
Fgr_new=zeros(sizeIm);
Ff_new=zeros(sizeIm);
Im_new=zeros(sizeIm); 
%Fgr(2:end-1,2:end-1,2:end-1)=((Im(3:end,2:end-1,2:end-1)-Im(1:end-2,2:end-1,2:end-1)).^2+(Im(2:end-1,3:end,2:end-1)-Im(2:end-1,1:end-2,2:end-1)).^2+(Im(2:end-1,2:end-1,3:end)-Im(2:end-1,2:end-1,1:end-2)).^2)./4;
Fgr(2:end-1,2:end-1,2:end-1)=((Im(3:end,2:end-1,2:end-1)-Im(2:end-1,2:end-1,2:end-1)).^2+(Im(1:end-2,2:end-1,2:end-1)-Im(2:end-1,2:end-1,2:end-1)).^2+...
    (Im(2:end-1,3:end,2:end-1)-Im(2:end-1,2:end-1,2:end-1)).^2+(Im(2:end-1,1:end-2,2:end-1)-Im(2:end-1,2:end-1,2:end-1)).^2+...
    (Im(2:end-1,2:end-1,3:end)-Im(2:end-1,2:end-1,2:end-1)).^2+(Im(2:end-1,2:end-1,1:end-2)-Im(2:end-1,2:end-1,2:end-1)).^2)./2;
Ff(Im(:)<=h)=alpha.*Im(Im(:)<=h)./h;
Ff(Im(:)>h)=alpha.*(1-Im(Im(:)>h))./(1-h);

figure, hold on
plot(0,mean(Fgr(:)),'b*')
plot(0,mean(Ff(:)),'r*')
plot(0,mean(Im(:)),'g*')
drawnow
legend('Fgr','Ff','Im')
    
for i=1:100
    Im_new(2:end-1,2:end-1,2:end-1)=Im(2:end-1,2:end-1,2:end-1)+((rand(sizeIm-2)>0.5).*2-1).*step_size;
    Im_new(Im_new(:)<0)=0;
    Im_new(Im_new(:)>1)=1;
    
    %Fgr_new(2:end-1,2:end-1,2:end-1)=((Im_new(3:end,2:end-1,2:end-1)-Im_new(1:end-2,2:end-1,2:end-1)).^2+(Im_new(2:end-1,3:end,2:end-1)-Im_new(2:end-1,1:end-2,2:end-1)).^2+(Im_new(2:end-1,2:end-1,3:end)-Im_new(2:end-1,2:end-1,1:end-2)).^2)./4;
    Fgr_new(2:end-1,2:end-1,2:end-1)=((Im(3:end,2:end-1,2:end-1)-Im_new(2:end-1,2:end-1,2:end-1)).^2+(Im(1:end-2,2:end-1,2:end-1)-Im_new(2:end-1,2:end-1,2:end-1)).^2+...
    (Im(2:end-1,3:end,2:end-1)-Im_new(2:end-1,2:end-1,2:end-1)).^2+(Im(2:end-1,1:end-2,2:end-1)-Im_new(2:end-1,2:end-1,2:end-1)).^2+...
    (Im(2:end-1,2:end-1,3:end)-Im_new(2:end-1,2:end-1,2:end-1)).^2+(Im(2:end-1,2:end-1,1:end-2)-Im_new(2:end-1,2:end-1,2:end-1)).^2)./2;
    Ff_new(Im_new(:)<=h)=alpha.*Im_new(Im_new(:)<=h)./h;
    Ff_new(Im_new(:)>h)=alpha.*(1-Im_new(Im_new(:)>h))./(1-h);
    
    ind=(rand(sizeIm)<exp(-(Fgr_new+Ff_new-Fgr-Ff)./T));
    ind(ind_fixed)=0;
    Im(ind)=Im_new(ind);
    Fgr(ind)=Fgr_new(ind);
    Ff(ind)=Ff_new(ind);
    
    plot(i,mean(Fgr(:)),'b*')
    plot(i,mean(Ff(:)),'r*')
    plot(i,mean(Im(:)),'g*')
    drawnow
end

Im=feval(cl,Im.*maxIm);



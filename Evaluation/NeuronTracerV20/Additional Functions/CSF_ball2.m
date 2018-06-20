% This is a multi-scale CSF - ball with a shell. Rball and Rshell can be vectors
% The volume of the ball is 1, and the volume of the shell is -1 

function Im_max = CSF_ball2(Im,Rball)

disp('CSF_ball filtering started.')
Rshell=(2^(1/3)-1).*Rball;
Rshell(Rshell<1)=1;
cl=class(Im);

if isempty(Rball) || prod(Rball)==0
    Im_max=Im;
else
    Im_max=zeros(size(Im));
    %Im_max=ones(size(Im)).*255;

    Im=double(Im);
    
    for k=1:length(Rball), [Rball(k),Rshell(k)]
        W=round(Rball(k)+Rshell(k))*2+1;
        center=(W+1)/2;
        CSF=zeros(W,W,W);
        [xx,yy,zz]=ind2sub_AS(size(CSF),[1:numel(CSF)]);
        r=((xx-center).^2+(yy-center).^2+(zz-center).^2).^0.5;
        CSF(r<=Rball(k))=1;
        CSF(r>Rball(k) & r<=(Rball(k)+Rshell(k)))=-1;
        CSF(r<=Rball(k))=CSF(r<=Rball(k)).*(1-sum(CSF(:))/sum(CSF(:)>0));
        CSF=CSF./(nnz(CSF(:)))^0.5;
        
        Im_temp = imfilter(Im,CSF,'replicate');
        %Im_temp=Im_temp./max(Im_temp(:)).*double(intmax(cl));
        
        Im_max=max(Im_max,Im_temp);
        %Im_max=min(Im_max,Im_temp);
    end
    Im_max=feval(cl,Im_max./max(Im_max(:)).*double(intmax(cl)));
end

disp('CSF_ball filtering is complete.')

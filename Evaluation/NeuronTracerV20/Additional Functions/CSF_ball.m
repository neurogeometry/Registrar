% This is a multi-scale CSF - ball with a shell. Rball and Rshell can be vectors
% The volume of the ball is 1, and the volume of the shell is -1 

function Im_max = CSF_ball(Im,Rball,Rshell)

disp('CSF_ball filtering started.')
cl=class(Im);

if isempty(Rball) || prod(Rball)==0 || length(Rball)~=length(Rshell)
    Im_max=Im;
else
    Im_max=zeros(size(Im));
    %Im_max=feval(cl,zeros(size(Im)));
    Im=double(Im);
    
    for k=1:length(Rball), [Rball(k),Rshell(k)]
        W=round(Rball(k)+Rshell(k))*2+1;
        center=(W+1)/2;
        CSF=zeros(W,W,W);
        [xx,yy,zz]=ind2sub_AS(size(CSF),[1:numel(CSF)]);
        r=((xx-center).^2+(yy-center).^2+(zz-center).^2).^0.5;
        CSF(r<=Rball(k))=1/nnz(r<=Rball(k)); %3/4/pi/Rball(k)^3;
        CSF(r>Rball(k) & r<=(Rball(k)+Rshell(k)))=-1/nnz(r>Rball(k) & r<=(Rball(k)+Rshell(k))); %-3/4/pi/((Rball(k)+Rshell(k))^3-Rball(k)^3);
        
        Im_temp = imfilter(Im,CSF,'replicate');
        %Im_temp=feval(cl,Im_temp./max(Im_temp(:)).*double(intmax(cl)));
        
        Im_max=max(Im_max,Im_temp);
    end
    Im_max=feval(cl,Im_max./max(Im_max(:)).*double(intmax(cl)));
end
disp('CSF_ball filtering is complete.')
























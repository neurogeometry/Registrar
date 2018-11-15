% This is a multi-scale G*cos CSF filter. S can be a vector 

function Im_max = CSF_cos(Im,S)

disp('CSF_cos filtering started.')
cl=class(Im);

if isempty(S) || prod(S)==0
    Im_max=Im;
else
    Im_max=zeros(size(Im));
    %Im_max=feval(cl,zeros(size(Im)));
    Im=double(Im);
    for k=1:length(S)
        W=round(3*S(k))*2+1;
        center=(W+1)/2;
        CSF=zeros(W,W,W);
        
        [xx,yy,zz]=ind2sub_AS(size(CSF),[1:numel(CSF)]);
        r=((xx-center).^2+(yy-center).^2+(zz-center).^2).^0.5;
        CSF(:)=exp(-r.^2./(2*S(k)^2)).*cos(r./S(k)); 
        CSF=CSF./sum(abs(CSF(:)));
        
        Im_temp = imfilter(Im,CSF,'replicate');
        %Im_temp=feval(cl,Im_temp./max(Im_temp(:)).*double(intmax(cl)));
        
        Im_max=max(Im_max,Im_temp);
        display(['scale = ', num2str(S(k))])
    end
    Im_max=feval(cl,Im_max./max(Im_max(:)).*double(intmax(cl)));
end
disp('CSF_cos filtering is complete.')
























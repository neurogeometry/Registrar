% This is a multi-scale CSF filter based on the Laplacian of G. S can be a vector 

function Im_max = CSF_d2G(Im,S)

disp('CSF_d2G filtering started.')
cl=class(Im);

if isempty(S) || prod(S)==0
    Im_max=Im;
else
    Im_max=zeros(size(Im));
    %Im_max=feval(cl,zeros(size(Im)));
    Im=double(Im);
    for k=1:length(S),S(k)
        W=round(3*S(k))*2+1;
        center=(W+1)/2;
        CSF=zeros(W,W,W);
        
        [xx,yy,zz]=ind2sub_AS(size(CSF),[1:numel(CSF)]);
        r=((xx-center).^2+(yy-center).^2+(zz-center).^2).^0.5;
        CSF(:)=exp(-r.^2./(2*S(k)^2)).*(1-(r./S(k)).^2./3); 
        CSF=CSF./sum(abs(CSF(:)));
        
        Im_temp = imfilter(Im,CSF,'replicate');
        %Im_temp=feval(cl,Im_temp./max(Im_temp(:)).*double(intmax(cl)));
        
        Im_max=max(Im_max,Im_temp);
    end
    Im_max=feval(cl,Im_max./max(Im_max(:)).*double(intmax(cl)));
end

disp('CSF_d2G filtering is complete.')
























% This filter eliminates image background, which is defined as (G*Im)<a*mean(Im). 
% R is defined as sqrt(2) times the SD of the Gaussian filter

function Im = Reduce_Background(Im,R,a)

disp('Background reduction started.')

if R>0 && a>=0
    cl=class(Im);
    
    sizeIm=size(Im);
    if length(sizeIm)==2
        sizeIm(3)=1;
    end

    W=ceil(2.5*R)*2+1;
    HW=(W-1)/2; % padding of Im
    center=(W+1)/2;
    %extra_padding=[0,0,0]; % extra padding to speed up fft
    extra_padding=[factorize(W-1+sizeIm(1)),factorize(W-1+sizeIm(2)),factorize(W-1+sizeIm(3))];
    
    fIm_pad=zeros(sizeIm+2.*HW+extra_padding,cl);
    fIm_pad(HW+1:HW+sizeIm(1),HW+1:HW+sizeIm(2),HW+1:HW+sizeIm(3))=Im;
    % replicate boundary conditions ???
    fIm_pad=fftn(fIm_pad);
    
    sizeIm_pad=size(fIm_pad);
    if length(sizeIm_pad)==2
        sizeIm_pad(3)=1;
    end
    %%%%%%%%%%%%%%%%% Gaussian Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    Kx=(2*pi*(0:sizeIm_pad(1)-1)./sizeIm_pad(1));
    Ky=(2*pi*(0:sizeIm_pad(2)-1)./sizeIm_pad(2));
    Kz=(2*pi*(0:sizeIm_pad(3)-1)./sizeIm_pad(3));
    Kxp=min(Kx,2*pi-Kx);
    Kyp=min(Ky,2*pi-Ky);
    Kzp=min(Kz,2*pi-Kz);
        
    fGx_pad=exp(-Kxp.^2.*R^2./4).*exp(1i.*Kx.*(center-1)); % conjugated fCSF
    fGy_pad=exp(-Kyp.^2.*R^2./4).*exp(1i.*Ky.*(center-1));
    fGz_pad=exp(-Kzp.^2.*R^2./4).*exp(1i.*Kz.*(center-1));
        
    Im_temp=real(ifftn(fIm_pad.*...
        repmat(transpose(fGx_pad),[1,sizeIm_pad(2),sizeIm_pad(3)]).*...
        repmat(fGy_pad,[sizeIm_pad(1),1,sizeIm_pad(3)]).*...
        repmat(reshape(fGz_pad,[1,1,sizeIm_pad(3)]),[sizeIm_pad(1),sizeIm_pad(2),1])));

%     if ~isa(Im,'double')
%         Im_temp=feval(cl,Im_temp./max(Im_temp(:)).*double(intmax(cl)));
%     end
    Im_temp=Im_temp(1:sizeIm(1),1:sizeIm(2),1:sizeIm(3));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Im=Im-feval(cl,a.*Im_temp);
    %Im(Im_temp<a)=0;
    %Im_temp=feval(cl,a.*Im_temp);
    %Im(Im<Im_temp)=0; %Im_temp(Im<Im_temp); 
end
disp('Background reduction is complete.')


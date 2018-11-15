% This is a Gaussian filter. 

function Im = CSF_G(Im,SD)

disp('CSF_G filter started.')
sizeIm=size(Im);

if SD>0
    cl=class(Im);
    Im=double(Im);

    W=ceil(3*SD)*2+1;
    HW=(W-1)/2; % padding of Im
    center=(W+1)/2;
    %extra_padding=[0,0,0]; % extra padding to speed up fft
    extra_padding=[factorize(W-1+size(Im,1)),factorize(W-1+size(Im,2)),factorize(W-1+size(Im,3))];
        
    Im_pad=zeros(size(Im)+2.*HW+extra_padding);
    Im_pad(HW+1:HW+size(Im,1),HW+1:HW+size(Im,2),HW+1:HW+size(Im,3))=Im;
    fIm_pad=fftn(Im_pad); % boundary conditions are not replicated ???
    
    %%%%%%%%%%%%%%%%% Gaussian Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    Kx=(2*pi*(0:size(Im_pad,1)-1)./size(Im_pad,1));
    Ky=(2*pi*(0:size(Im_pad,2)-1)./size(Im_pad,2));
    Kz=(2*pi*(0:size(Im_pad,3)-1)./size(Im_pad,3));
    Kxp=min(Kx,2*pi-Kx);
    Kyp=min(Ky,2*pi-Ky);
    Kzp=min(Kz,2*pi-Kz);
        
    fGx_pad=exp(-Kxp.^2.*SD^2./2).*exp(1i.*Kx.*(center-1)); % conjugated fCSF
    fGy_pad=exp(-Kyp.^2.*SD^2./2).*exp(1i.*Ky.*(center-1));
    fGz_pad=exp(-Kzp.^2.*SD^2./2).*exp(1i.*Kz.*(center-1));
        
    fG_pad=repmat(transpose(fGx_pad),[1,size(Im_pad,2),size(Im_pad,3)]).*...
        repmat(fGy_pad,[size(Im_pad,1),1,size(Im_pad,3)]).*...
        repmat(reshape(fGz_pad,[1,1,size(Im_pad,3)]),[size(Im_pad,1),size(Im_pad,2),1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    Im=ifftn(fIm_pad.*fG_pad);
    Im=real(Im(1:sizeIm(1),1:sizeIm(2),1:sizeIm(3)));
    Im=feval(cl,Im./max(Im(:)).*double(intmax(cl)));    
end
disp('CSF_G filtering is complete.')




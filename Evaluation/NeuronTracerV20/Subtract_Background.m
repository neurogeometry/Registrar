% This filter subtracts image background, which is defined as a*(G*Im).
% R is defined as sqrt(2) times the SD of the Gaussian filter, G

function Im = Subtract_Background(Im,R,a)

disp('Background subtraction started.')

if R>0 && a>0
    cl=class(Im);
    
%     b=1;
%     if isa(Im,'uint8')
%         b=1;
%     elseif isa(Im,'uint16')
%         b=2;
%     elseif isa(Im,'uint32')
%         b=4;
%     elseif isa(Im,'uint64')
%         b=8;
%     elseif isa(Im,'double')
%         b=8;
%     end
    
    sizeIm=size(Im);
    if length(sizeIm)==2
        sizeIm(3)=1;
    end
    
    W=ceil(2.5*R)*2+1;
    HW=(W-1)/2; % padding of Im
    center=(W+1)/2;
    %extra_padding=[0,0,0]; % extra padding to speed up fft
    extra_padding=[factorize(W-1+sizeIm(1)),factorize(W-1+sizeIm(2)),factorize(W-1+sizeIm(3))];
    
    out = imaqmex_1('memorystatus');
    %out = imaqmem;
    mem_left = out.AvailPhys;
    Byte_needed=(8*2*2*1.5+1)*prod(sizeIm+2.*HW+extra_padding);
    
    if(Byte_needed > mem_left)
        disp('Insufficient amount of available memory to run this plugin. Consider reducing the image.');
    else
        %fIm_pad=zeros(sizeIm+2.*HW+extra_padding,cl);
        %fIm_pad(HW+1:HW+sizeIm(1),HW+1:HW+sizeIm(2),HW+1:HW+sizeIm(3))=Im;
        
        fIm_pad=padarray(Im,[HW HW HW],'symmetric'); %'replicate'
        fIm_pad=padarray(fIm_pad,extra_padding,'post');
        
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
        
        Im_temp=Im_temp(1:sizeIm(1),1:sizeIm(2),1:sizeIm(3));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Im=Im-feval(cl,a.*Im_temp);
    end
end
disp('Background subtraction is complete.')


% This is a fast multi-scale LoG filter. R is the scale,
% defined as sqrt(2) times the SD of the best LoG filter

function Im_max = Multi_Scale_LoG_2D(Im,R_min,R_step,R_max)

disp('Multi-scale LoG filtering started.')
cl=class(Im);

% b=1;
% if isa(Im,'uint8')
%     b=1;
% elseif isa(Im,'uint16')
%     b=2;
% elseif isa(Im,'uint32')
%     b=4;
% elseif isa(Im,'uint64')
%     b=8;
% elseif isa(Im,'double')
%     b=8;
% end

sizeIm=size(Im);
if length(sizeIm)==2
    sizeIm(3)=1;
end

if R_min<0.1
    R_min=0.1;
end

s=R_min:R_step:R_max;

if isempty(s)
    Im_max=Im;
else
    del=0.01*R_min;
    
    s0=max(s);
    W=ceil(2.5*s0)*2+1;
    HW=(W-1)/2; % padding of Im
    center=(W+1)/2;
    %extra_padding=[0,0]; % extra padding to speed up fft
    extra_padding=[factorize(W-1+sizeIm(1)),factorize(W-1+sizeIm(2))];
    
    out = imaqmex_1('memorystatus');
    %out = imaqmem;
    mem_left = out.AvailPhys;
    Byte_needed=(8*2*2*2+1)*prod(sizeIm(1:2)+2.*HW+extra_padding);
    
    if(Byte_needed > mem_left)
        disp('Insufficient amount of available memory to run this plugin. Consider reducing the image.');
        Im_max=Im; 
    else
        Im_max=zeros(sizeIm,cl);
        
        %fIm_pad=zeros(sizeIm(1:2)+2.*HW+extra_padding,cl);
        %fIm_pad(HW+1:HW+sizeIm(1),HW+1:HW+sizeIm(2))=Im;
        
        fIm_pad=padarray(Im,[HW HW],'symmetric'); %'replicate'
        fIm_pad=padarray(fIm_pad,extra_padding,'symmetric','post');
        
        % replicate boundary conditions ???
        clear Im
        fIm_pad=fft2(fIm_pad);
        sizeIm_pad=size(fIm_pad);
        
        Kx=(2*pi*(0:sizeIm_pad(1)-1)./sizeIm_pad(1));
        Ky=(2*pi*(0:sizeIm_pad(2)-1)./sizeIm_pad(2));
        Kxp=min(Kx,2*pi-Kx);
        Kyp=min(Ky,2*pi-Ky);
        
        for k=1:length(s)
            %%%%%%%%%%%%%%%%% LoG Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fCSFx_pad1=exp(-Kxp.^2.*(s(k)-del)^2./4).*exp(1i.*Kx.*(center-1)); % conjugated fCSF
            fCSFy_pad1=exp(-Kyp.^2.*(s(k)-del)^2./4).*exp(1i.*Ky.*(center-1));
            
            fCSFx_pad2=exp(-Kxp.^2.*(s(k)+del)^2./4).*exp(1i.*Kx.*(center-1)); % conjugated fCSF
            fCSFy_pad2=exp(-Kyp.^2.*(s(k)+del)^2./4).*exp(1i.*Ky.*(center-1));
            
            Im_temp=real(ifft2(fIm_pad.*(repmat(transpose(fCSFx_pad1),[1,sizeIm_pad(2)]).*repmat(fCSFy_pad1,[sizeIm_pad(1),1])-...
                repmat(transpose(fCSFx_pad2),[1,sizeIm_pad(2)]).*repmat(fCSFy_pad2,[sizeIm_pad(1),1]))./(2*del/s(k)))); % correct up to a numerical factor
            
            Im_temp(Im_temp<0)=0;
            if ~isa(Im_max,'double')
                Im_temp=feval(cl,Im_temp./max(Im_temp(:)).*double(intmax(cl)));
            end
            Im_temp=Im_temp(1:sizeIm(1),1:sizeIm(2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % combining different scales
            if k==1
                Im_max=Im_temp;
                clear Im_temp
            else
                ind=(Im_temp > Im_max) & (Im_temp > Simple_Filters(Im_max,'Max',[(2*s(k-1)+1).*ones(1,2),1]));
                Im_max(ind)=Im_temp(ind);
                clear Im_temp ind
            end
            
            display(['scale = ', num2str(s(k))])
        end
    end
end

disp('Multi-scale LoG filtering is complete.')




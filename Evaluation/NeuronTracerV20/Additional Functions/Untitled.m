S=2.3;
W=ceil(3*S)*2+1;
HW=(W-1)/2;
center=(W+1)/2;
CSF=[exp(-(-HW:HW).^2./2./S^2)./(2*pi*S^2)^0.5];
% N=length(CSF);
% K=(2*pi*(0:N-1)./N);
% Kp=min(K,2*pi-K);
% fCSF=exp(-Kp.^2.*S^2./2-1i.*K.*(center-1)); %works even if CSF is padded from the right
%plot(real(fft(CSF))), hold on, plot(real(fCSF),'r')

x=1:100;
Im=rand(1,length(x)); %1./x+1;
extra_pad=0;
pad=HW+extra_pad;  
Im_pad=zeros(1,length(Im)+2*pad);
Im_pad(pad+1:end-pad)=Im; 
%Im_pad(1:pad)=Im(pad:-1:1); Im_pad(length(Im_pad)-pad+1:length(Im_pad))=Im(length(Im):-1:length(Im)-pad+1);

CSF_pad=zeros(1,length(Im_pad));
CSF_pad(1:length(CSF))=CSF;
N=length(CSF_pad);
K=(2*pi*(0:N-1)./N);
Kp=min(K,2*pi-K);
fCSF_pad=exp(-Kp.^2.*S^2./2-1i.*K.*(center-1)); %works even if CSF is padded from the right


%Im_rep=imfilter(Im,CSF,'replicate');
%Im_ref=imfilter(Im,CSF,'symmetric');
Im_0=imfilter(Im,CSF);
Imp=ifft(fft(Im_pad).*conj(fCSF_pad));
Imp=real(Imp(extra_pad+1:extra_pad+length(Im)));

plot(x,Im)
hold on
plot(x,Imp,'k-')
%plot(x,Im_rep,'r-')
%plot(x,Im_ref,'g-')
plot(x,Im_0,'m-')



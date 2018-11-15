
% FIJI Registration -> Elastic -> Elastic stack alignement
clear all;
rootpath = 'E:\SliceRegistration\DL083B001G_Mess_Elastic_Fiji\';

k = 1;
for i =0:97
    
    IM1 = imread([rootpath,'elastic-000',num2str(i,'%02.f'),'.tif']);
    IM2 = imread([rootpath,'elastic-000',num2str(i+2,'%02.f'),'.tif']);
    
%     imshowpair(uint16(IM1).*600,uint16(IM2).*600,'Scaling','independent')
    cor3 = corrcoef(double(IM1),double(IM2));
    if max(max(isnan(cor3))) == 0
        correl(k) = cor3(1,2)
        k = k + 1;
    end
end
mean(correl)
std(correl)









% FIJI SIFT RIGID
clear all;
rootpath = 'E:\SliceRegistration\DL083B001G_Mess_SIFT_Fiji\Affine\'; % RIGID % Translation

k = 1;
for i =0:98
    
    IM1 = imread([rootpath,'Aligned 101 of 10100',num2str(i,'%02.f'),'.tif']);
    IM2 = imread([rootpath,'Aligned 101 of 10100',num2str(i+1,'%02.f'),'.tif']);
    
    % imshowpair(uint16(IM1).*600,uint16(IM2).*600,'Scaling','independent')
    cor3 = corrcoef(double(IM1),double(IM2));
    if max(max(isnan(cor3))) == 0
        correl(k) = cor3(1,2)
        k = k + 1;
    end
end
mean(correl)
std(correl)














% IM1 = imread('E:\SliceRegistration\DL083B001G_Mess_bUnwarpJ\Registered Target Image0000.tif');
% IM2 = imread('E:\SliceRegistration\DL083B001G_Mess_bUnwarpJ\Registered Target Image0001.tif');

% Fiji bUnwarp
clear all;
rootpath = 'E:\SliceRegistration\DL083B001G_Mess_bUnwarpJ_1-2';%     E:\SliceRegistration\DL083B001G_Mess_bUnwarpJ_Non_Rigid


for i =0:99
    
IM1 = imread([rootpath,'\',num2str(i),'\Registered Target Image0000.tif']);
IM2 = imread([rootpath,'\',num2str(i),'\Registered Target Image0001.tif']);
% imshowpair(uint16(IM1).*600,uint16(IM2).*600,'Scaling','independent')
cor3 = corrcoef(double(IM1),double(IM2));
correl(i+1) = cor3(1,2)
end
mean(correl)
std(correl)

figure,
boxplot([correl(:)],'Whisker',inf)
    axis square, box on

% correl(12)
% for i=16:99
%     
%     
%     mkdir([rootpath,'\',num2str(i)])
%     
% end
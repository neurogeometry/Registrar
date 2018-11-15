function im = appendimages(image1, image2,Direction)
if strcmp(Direction,'horizontal')
    rows1 = size(image1,1);
    rows2 = size(image2,1);
    
    if (rows1 < rows2)
        image1(rows2,1) = 0;
    else
        image2(rows1,1) = 0;
    end
%     padd = ones(size(image1,2),45);
%     padd(:,:) = max(image1(:));
%     im = [image1 padd image2];
     im = [image1 image2];
else
    cols1 = size(image1,2);
    cols2 = size(image2,2);
    
    if (cols1 < cols2)
        image1(cols1,2) = 0;
    else
        image2(cols2,2) = 0;
    end
    
    im = [image1; image2];
end





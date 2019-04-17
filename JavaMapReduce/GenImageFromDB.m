SaveFolder = 'C:\Users\Seyed\Documents\DatasetTests\NCTWeb\NCTWeb\NCT-Web\src\main\webapp\resources\core\ImageMouseLight\';
StackClass='uint8';
z_level = 1;
imageID = 3;
 conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:','E:/nctracer_Joe_AllTiles_AllZoomLevels.db']);
%conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:','F:/nctracer_MouseLight_200_200_100.db']);
% x = conn.Handle;

SqlStr = ['SELECT x,y,z FROM pix where image_id = ',num2str(imageID),' and z_level = ',num2str(z_level)];
result = fetch(conn,SqlStr);


for i = 1:size(result,1)
    x = result{i,1};
    y = result{i,2};
    z = result{i,3};
    SqlStr = ['SELECT pixels,x_dim,y_dim,z_dim FROM pix where x = ',num2str(x),' and y = ',num2str(y),' and z = ',num2str(z),' and z_level = ',num2str(z_level),' and image_id = ',num2str(imageID)];
    result1 = fetch(conn,SqlStr);
    Tile = reshape(typecast(result1{1},StackClass),result1{2},result1{3},result1{4});
    Tile = imrotate(Tile,90);
    TileName = [num2str(x),'_',num2str(y),'_',num2str(z)];
    mkdir([SaveFolder,TileName]);
    for zz = 0:size(Tile,3)-1       
        TileName_1 = [sprintf('%04d', zz),'.png'];
        imwrite(Tile(:,:,zz+1),[SaveFolder,TileName,'/',TileName_1],'png','Compression','jpeg');
    end
end
close(conn)
% figure,imshow(max(IM,[],3),[0 max(IM(:))])
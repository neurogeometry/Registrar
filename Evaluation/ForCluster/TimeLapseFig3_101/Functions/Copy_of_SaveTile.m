function SaveTile(usedDB,Tile_new,image_id,x,TilePositions,TilePositions_new,TileSize_org,i,t,z_level,N_tiles,SaveFolder,N_tiles_new,namecount)
Tile_glancer = [];
TileName=num2str(namecount,['%0',num2str(fix(log10(prod(N_tiles)*prod(N_tiles_new)))+1),'.0f']);
if usedDB == 1
    raw_im = reshape(Tile_new,1,[]);
    insertcommand = 'INSERT INTO pix (image_id,name,x,y,z,x_dim,y_dim,z_dim,pixels,zoom_out) values (?,?,?,?,?,?,?,?,?,?)';
    StatementObject = x.prepareStatement(insertcommand);
    StatementObject.setObject(1,image_id);
    StatementObject.setObject(2,TileName);
    StatementObject.setObject(3,TilePositions(i,1)+TilePositions_new(t,1)-1)
    StatementObject.setObject(4,TilePositions(i,2)+TilePositions_new(t,2)-1)
    StatementObject.setObject(5,TilePositions(i,3)+TilePositions_new(t,3)-1)
    StatementObject.setObject(6,TileSize_org(1));
    StatementObject.setObject(7,TileSize_org(2));
    StatementObject.setObject(8,TileSize_org(3));
    StatementObject.setObject(9,raw_im);
    StatementObject.setObject(10,z_level);
    StatementObject.execute;
    close(StatementObject);
    disp(['Tile ',num2str(i),'/',num2str(prod(N_tiles)),' inserted.']);
elseif usedDB == 2
    %                 TileName=num2str(i,['%0',num2str(fix(log10(prod(N_tiles)))+1),'.0f']);
    mkdir([SaveFolder,TileName]);
%     saveastiff(im2uint8(Tile_new), [SaveFolder,TileName,'\',TileName,'.tif']);
    options.overwrite = 1;
    saveastiff(Tile_new, [SaveFolder,TileName,'\',TileName,'.tif'],options);
    disp(['Tile ',num2str(i),'/',num2str(prod(N_tiles)),' created.']);
elseif usedDB == 3 % Neuroglancer
    %                 TileName=num2str(i,['%0',num2str(fix(log10(prod(N_tiles)))+1),'.0f']);
%     for ii = 1:size(Tile_new,3)
%         Tile_glancer = [Tile_glancer;Tile_new(:,:,ii)];
%     end
    C = permute(Tile_new,[1 3 2]);
    Tile_glancer = reshape(C,[],size(Tile_new,2),1);

    
    xstart = TilePositions(i,2)+TilePositions_new(t,2)-2;
    xend = xstart+size(Tile_new,1);
    ystart = TilePositions(i,1)+TilePositions_new(t,1)-2;
    yend = ystart+size(Tile_new,2);
    zstart = TilePositions(i,3)+TilePositions_new(t,3)-2;
    zend = zstart+size(Tile_new,3);
    
    TileName = [num2str(xstart),'-',num2str(xend),'_',num2str(ystart),'-',num2str(yend),'_',num2str(zstart),'-',num2str(zend)];
    %                 saveastiff(im2uint8(Tile_glancer), [SaveFolder,'image\',TileName]);
    imwrite(Tile_glancer,[SaveFolder,'image\',TileName],'jpg');
    Tile_glancer = [];
    disp(['Tile ',num2str(i),'/',num2str(prod(N_tiles)),' created.']);
elseif usedDB == 4 % Nifti
    xstart = TilePositions(i,2)+TilePositions_new(t,2)-2;
    ystart = TilePositions(i,1)+TilePositions_new(t,1)-2;
    zstart = TilePositions(i,3)+TilePositions_new(t,3)-2;
    mkdir([SaveFolder,'Nifti\',TileName]);
    Nifti_Save(Tile_new,xstart,ystart,zstart,TileName,[SaveFolder,TileName]);
    
elseif usedDB == 5 % HDF5 - Big Data Viewer (Fiji)
    
    resolutions =  ones(prod(N_tiles)*prod(N_tiles_new),1);
    subdivisions = repmat(16,prod(N_tiles)*prod(N_tiles_new),1);
    DbName = [SaveFolder,'HDf5\export1.h5'];
    if namecount < 11
        S =num2str(namecount-1,['%0',num2str(2),'.0f'])
    else
        S =num2str(namecount-1)
    end
    h5create(DbName,['/t00000/s',S,'/0/cells'],size(Tile_new),'ChunkSize',size(Tile_new));
    h5write(DbName, ['/t00000/s',S,'/0/cells'], im2uint8(Tile_new));
    fileattrib(DbName,'+w');
    h5writeatt(DbName,['/t00000/s',S,'/0/cells'],'element_size_um',[0.5,0.5,0.5]);
    
    h5create(DbName,['/s',S,'/resolutions'],size(resolutions),'ChunkSize',size(resolutions));
    h5write(DbName, ['/s',S,'/resolutions'], resolutions);
    fileattrib(DbName,'+w');
    h5writeatt(DbName,['/s',S,'/resolutions'],'element_size_um',[0.5,0.5,0.5]);
    
    h5create(DbName,['/s',S,'/subdivisions'],size(subdivisions),'ChunkSize',size(subdivisions));
    h5write(DbName, ['/s',S,'/subdivisions'], subdivisions);
    fileattrib(DbName,'+w');
    h5writeatt(DbName,['/s',S,'/subdivisions'],'element_size_um',[0.5,0.5,0.5]);
elseif usedDB == 6 % CATMAID - TrackEM
    xstart = (TilePositions(i,1)+TilePositions_new(t,1)-2)/TileSize_org(1);
    ystart = (TilePositions(i,2)+TilePositions_new(t,2)-2)/TileSize_org(2);
    zstart = TilePositions(i,3)+TilePositions_new(t,3)-2;
    
    
    for z = 0:size(Tile_new,3)-1
        folder = zstart+z;
        if ~exist([SaveFolder,'CATMAID\',num2str(folder)],'dir')
            mkdir([SaveFolder,'CATMAID\',num2str(folder)]);
        end
        TileName = [num2str(xstart),'_',num2str(ystart),'_0.jpg'];
        
        SubtileTile_CATMAID = Tile_new(:,:,z+1);
        imwrite(SubtileTile_CATMAID,[SaveFolder,'CATMAID\',num2str(folder),'\',TileName],'jpg');
        
        
    end
    
end
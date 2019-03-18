function SaveTile(usedDB,Tile,image_id,x,TilePositions,TileSize,z_level,SaveFolder,N_tiles)
TileName=[num2str(TilePositions(1)),'_',num2str(TilePositions(2)),'_',num2str(TilePositions(3))];
if usedDB == 1
    raw_im = reshape(Tile,1,[]);
    insertcommand = 'INSERT INTO pix (image_id,x,y,z,x_dim,y_dim,z_dim,pixels,zoom_out) values (?,?,?,?,?,?,?,?,?)';
    StatementObject = x.prepareStatement(insertcommand);
    StatementObject.setObject(1,image_id);
    StatementObject.setObject(2,TilePositions(1)-1)
    StatementObject.setObject(3,TilePositions(2)-1)
    StatementObject.setObject(4,TilePositions(3)-1)
    StatementObject.setObject(5,TileSize(1));
    StatementObject.setObject(6,TileSize(2));
    StatementObject.setObject(7,TileSize(3));
    StatementObject.setObject(8,raw_im);
    StatementObject.setObject(9,z_level);
    StatementObject.execute;
    close(StatementObject);
elseif usedDB == 2
    mkdir([SaveFolder,TileName]);
    options.overwrite = 1;
    saveastiff(Tile, [SaveFolder,TileName,'/',TileName,'.tif'],options);
    
    % for save as JPEG (Like Neuroglancer) for Joe
    C1 = permute(Tile,[1 3 2]);
    Tile_glancer1 = reshape(C1,[],size(Tile,2),1);
    xstart1 = TilePositions(2)-1;
    xend1 = xstart1+size(Tile,1);
    ystart1 = TilePositions(1)-1;
    yend1 = ystart1+size(Tile,2);
    zstart1 = TilePositions(3)-1;
    zend1 = zstart1+size(Tile,3);
    TileName1 = [num2str(xstart1),'-',num2str(xend1),'_',num2str(ystart1),'-',num2str(yend1),'_',num2str(zstart1),'-',num2str(zend1)];
    imwrite(Tile_glancer1,[SaveFolder,TileName1,'.jpg'],'jpg');
elseif usedDB == 3 % Neuroglancer
    C = permute(Tile,[1 3 2]);
    Tile_glancer = reshape(C,[],size(Tile,2),1);
    xstart = TilePositions(2)-1;
    xend = xstart+size(Tile,1);
    ystart = TilePositions(1)-1;
    yend = ystart+size(Tile,2);
    zstart = TilePositions(3)-1;
    zend = zstart+size(Tile,3);
    
    TileName = [num2str(xstart),'-',num2str(xend),'_',num2str(ystart),'-',num2str(yend),'_',num2str(zstart),'-',num2str(zend)];
    imwrite(Tile_glancer,[SaveFolder,'image/',TileName],'jpg');
elseif usedDB == 4 % Nifti
    xstart = TilePositions(2)-2;
    ystart = TilePositions(1)-2;
    zstart = TilePositions(3)-2;
    mkdir([SaveFolder,'Nifti/',TileName]);
    Nifti_Save(Tile,xstart,ystart,zstart,TileName,[SaveFolder,TileName]);
elseif usedDB == 5 % HDF5 - Big Data Viewer (Fiji)
    resolutions =  ones(prod(N_tiles),1);
    subdivisions = repmat(16,prod(N_tiles),1);
    DbName = [SaveFolder,'HDf5/export1.h5'];
    h5create(DbName,['/t00000/s',TileName,'/0/cells'],size(Tile),'ChunkSize',size(Tile));
    h5write(DbName, ['/t00000/s',TileName,'/0/cells'], im2uint8(Tile));
    fileattrib(DbName,'+w');
    h5writeatt(DbName,['/t00000/s',TileName,'/0/cells'],'element_size_um',[0.5,0.5,0.5]);
    
    h5create(DbName,['/s',TileName,'/resolutions'],size(resolutions),'ChunkSize',size(resolutions));
    h5write(DbName, ['/s',TileName,'/resolutions'], resolutions);
    fileattrib(DbName,'+w');
    h5writeatt(DbName,['/s',TileName,'/resolutions'],'element_size_um',[0.5,0.5,0.5]);
    
    h5create(DbName,['/s',TileName,'/subdivisions'],size(subdivisions),'ChunkSize',size(subdivisions));
    h5write(DbName, ['/s',TileName,'/subdivisions'], subdivisions);
    fileattrib(DbName,'+w');
    h5writeatt(DbName,['/s',TileName,'/subdivisions'],'element_size_um',[0.5,0.5,0.5]);
elseif usedDB == 6 % CATMAID - TrackEM
    xstart = (TilePositions(1)-2)/TileSize(1);
    ystart = (TilePositions(2)-2)/TileSize(2);
    zstart = TilePositions(3)-2;

    for z = 0:size(Tile,3)-1
        folder = zstart+z;
        if ~exist([SaveFolder,'CATMAID/',num2str(folder)],'dir')
            mkdir([SaveFolder,'CATMAID/',num2str(folder)]);
        end
        TileName = [num2str(xstart),'_',num2str(ystart),'_0.jpg'];
        
        SubtileTile_CATMAID = Tile(:,:,z+1);
        imwrite(SubtileTile_CATMAID,[SaveFolder,'CATMAID/',num2str(folder),'/',TileName],'jpg');
    end
elseif usedDB == 7
    mkdir([SaveFolder,'HDF5']);
    
    hdf5write([SaveFolder,'HDF5','temp_',TileName,'.h5'], '/dataset1', Tile);
end
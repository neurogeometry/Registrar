% This function creates seamless image tiles from the dataset of original image stacks and stack transformations.

function Retiling(StackList,StackPositions,StackSizes,T,TileSize,DataFolder,usedDB)
% addpath('../Functions');
parameters;
TileSize_org = TileSize;
% InfoImage=imfinfo(char(StackList(1,1)));
% M_tiles = ceil(InfoImage(1).Height/TileSize(1));
% TileSize = TileSize.*2;

image_id = 1;
x = [];
z_level = 1;

[filepath,~,ext] = fileparts(char(StackList(1,1)));

if params.RE.useHDF5
     X = hdf5read([DataFolder,'\tmp\temp_',num2str(1),'.h5'], '/dataset1');
     MaxIntensityValue = max(X(:));
    StackClass = class(hdf5read([DataFolder,'\tmp\temp_',num2str(1),'.h5'], '/dataset1'));  
else  
    if size(ext,2)>1
        InfoImage=imfinfo(char(StackList(1,1)));
        MaxIntensityValue = InfoImage(1).MaxSampleValue;
        StackClass = class(imread(char(StackList(1,1))));
    else
        allfiles = dir(filepath);
        InfoImage=imfinfo(char([allfiles(3).folder,'\',allfiles(3).name]));
        MaxIntensityValue = InfoImage(1).MaxSampleValue;
        StackClass = class(imread(char([allfiles(3).folder,'\',allfiles(3).name])));
    end
end

% StackClass='uint16';
% TileSize_org = TileSize;
% TileSize = StackSizes(1,:).*2;
SpecimenName = extractBefore(extractAfter(DataFolder,"MicroscopeFiles\Results-"),'_StackList');
if usedDB == 1
    %     DBFile = [pwd,'\',DataFolder,'\nctracer.db'];
    DBFile = 'E:\TilesCreation\NCTracerWeb\New\NCtracerWeb-master\NCtracerWeb-master\NCT-Web\data\db\nctracer.db';
    % for connection help: https://www.mathworks.com/help/database/ug/sqlite-jdbc-windows.html#bt8kopj-1
    if exist(DBFile) > 0
        conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:',DBFile]);
        conn.Message
        x = conn.Handle;     
    else
        
        conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:',DBFile]);
        conn.Message
        x = conn.Handle;
        sqlquery_createEDGE = 'CREATE TABLE edge (id integer primary key autoincrement not null, trace_id integer not null, v1_id integer not null, v2_id integer not null)';
        curs = exec(conn,sqlquery_createEDGE);
        close(curs);
        
        sqlquery_createIMAGE = 'CREATE TABLE image (id integer primary key autoincrement not null, name text not null)';
        curs = exec(conn,sqlquery_createIMAGE);
        close(curs);
        
        sqlquery_createPIX = 'CREATE TABLE "pix" ( `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT, `name` TEXT NOT NULL,`image_id` integer NOT NULL, `x` integer NOT NULL, `y` integer NOT NULL, `z` integer NOT NULL, `x_dim` integer NOT NULL, `y_dim` integer NOT NULL, `z_dim` integer NOT NULL, `pixels` blob NOT NULL,`zoom_out` INTEGER NOT NULL)';
        curs = exec(conn,sqlquery_createPIX);
        close(curs);
        
        sqlquery_createTRACE = 'CREATE TABLE trace (id integer primary key autoincrement not null, name text not null, image_id integer not null)';
        curs = exec(conn,sqlquery_createTRACE);
        close(curs);
        
        sqlquery_createVERTEX = 'CREATE TABLE vertex (id integer primary key autoincrement not null, trace_id integer not null, x integer not null, y integer not null, z integer not null)';
        curs = exec(conn,sqlquery_createVERTEX);
        close(curs);
        
        
        sqlquery_createIndex = 'create index ix_image_id on image (id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_pix_id on pix (id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_pix_image_id on pix (image_id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_pix_x on pix (x)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_pix_y on pix (y)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_pix_z on pix (z)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_pix_zoom on pix (zoom_out)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        
        sqlquery_createIndex = 'create index ix_trace_id on trace (id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_vertex_id on vertex (id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_vertex_trace_id on vertex (trace_id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_edge_id on edge (id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_edge_trace_id on edge (trace_id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_edge_v1_id on edge (v1_id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);
        
        sqlquery_createIndex = 'create index ix_edge_v2_id on edge (v2_id)';
        curs = exec(conn,sqlquery_createIndex);
        close(curs);

    end
    
    insertcommand_IMAGE = 'INSERT INTO image (name) values (?)';
    StatementObject = x.prepareStatement(insertcommand_IMAGE);
    StatementObject.setObject(1,SpecimenName);
    StatementObject.execute
    close(StatementObject);
    
    q2 = ("select seq from sqlite_sequence where name='image'");
    curs = exec(conn, q2);
    curs = fetch(curs);
    image_id = curs.Data{1}
    z_level = 1;
end

SaveFolder=[DataFolder,params.RE.savefolder,'Zoom1\'];
% for test
SaveFolder = 'E:\TilesCreation\Data\MouselightFull_Neuroglancer\';
N_stacks=length(StackList);
if usedDB == 3
    mkdir([SaveFolder,'image']);
elseif usedDB == 4
    mkdir([SaveFolder,'Nifti']);
elseif usedDB == 5
    mkdir([SaveFolder,'HDF5']); 
elseif usedDB == 6
    mkdir([SaveFolder,'CATMAID']); 
end
StackHW=StackSizes./2;

if strcmp(T.transform,'Translation')
    % Transform the binding boxes of all stacks
    b=T.b';
    
    Verts=zeros(8,3,N_stacks);
    Verts_transformed=zeros(8,3,N_stacks);
    for i=1:N_stacks
        Verts(:,:,i)=ones(8,1)*StackPositions(i,:)-1+[1,1,1;StackSizes(i,1),1,1;1,StackSizes(i,2),1;StackSizes(i,1),StackSizes(i,2),1;...
            1,1,StackSizes(i,3);StackSizes(i,1),1,StackSizes(i,3);1,StackSizes(i,2),StackSizes(i,3);StackSizes(i,1),StackSizes(i,2),StackSizes(i,3)];
        Verts_transformed(:,:,i)=Verts(:,:,i)+b(i,:);
    end
    
    % Find the extent of the transformed brain space
    Max=max(max(Verts_transformed,[],3),[],1);
    Min=min(min(Verts_transformed,[],3),[],1);
    Global_shift=1-Min;
    
    StackPositionsTransformed=round(StackPositions+b+ones(N_stacks,1)*Global_shift);
    StackCentersTransformed=StackPositionsTransformed+(StackSizes-1)./2;
    
    % Tile this space
    N_tiles=ceil(ceil(Max-Min)./TileSize);
    [xx,yy,zz]=ndgrid(0:N_tiles(1)-1,0:N_tiles(2)-1,0:N_tiles(3)-1);
    TilePositions=1+[xx(:).*TileSize(1),yy(:).*TileSize(2),zz(:).*TileSize(3)];
    TileHW=TileSize./2;
    namecount = 1;
    for i=1:prod(N_tiles)
        % For a given tile, find stacks that overlap it after translation
        TileCenter=TilePositions(i,:)+(TileSize-1)./2;
        StackInd=find(abs(StackCentersTransformed(:,1)-TileCenter(1))<StackHW(:,1)+TileHW(1) & abs(StackCentersTransformed(:,2)-TileCenter(2))<StackHW(:,2)+TileHW(2) & abs(StackCentersTransformed(:,3)-TileCenter(3))<StackHW(:,3)+TileHW(3));
        Tile=zeros(TileSize,StackClass);
        
        for j=1:length(StackInd)
            % Write overlapping stack intensity into the tile

            if params.RE.useHDF5
                X = hdf5read([DataFolder,'\tmp\temp_',num2str(StackInd(j)),'.h5'], '/dataset1');
            else
                pth=StackList{StackInd(j)}(1:find(StackList{StackInd(j)}=='\',1,'last'));
                file_list=StackList{StackInd(j)}(find(StackList{StackInd(j)}=='\',1,'last')+1:end);
               % X=im2uint8(ImportStack([pth,file_list]));
                             X=ImportStack([pth,file_list]);
            end          
%                          figure();imshow(max(X,[],3),[0 max(X(:))]);
            %             [X,~,~]=ImportStackJ(pth,{file_list},1,1,1,'Max');
            
            TileStart=max(StackPositionsTransformed(StackInd(j),:)-TilePositions(i,:)+1,1);
            TileEnd=min(StackPositionsTransformed(StackInd(j),:)+StackSizes(StackInd(j),:)-TilePositions(i,:),TileSize);
            StackStart=max(1,TilePositions(i,:)-StackPositionsTransformed(StackInd(j),:)+1);
            StackEnd=min(StackSizes(StackInd(j),:),TilePositions(i,:)+TileSize-StackPositionsTransformed(StackInd(j),:));
            Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=...
                max(Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3)),...
                X(StackStart(1):StackEnd(1),StackStart(2):StackEnd(2),StackStart(3):StackEnd(3)));
%             figure();imshow(max(Tile,[],3),[0 max(Tile(:))]);

        end

        Tile = uint8(double(Tile)./double(MaxIntensityValue)*255);
%         figure();imshow(max(Tile,[],3),[0 max(Tile(:))]);
 
        N_tiles1 = 1;
        Verts=zeros(8,3,N_tiles1);
        
        Verts(:,:)=ones(8,1)*TilePositions(i,:)-1+[1,1,1;TileSize(1),1,1;1,TileSize(2),1;TileSize(1),TileSize(2),1;...
            1,1,TileSize(3);TileSize(1),1,TileSize(3);1,TileSize(2),TileSize(3);TileSize(1),TileSize(2),TileSize(3)];
        
        % Find the extent of the transformed brain space
        Max=max(max(Verts,[],3),[],1);
        Min=min(min(Verts,[],3),[],1);
        
        % Tile this space
        N_tiles_new=ceil(ceil(Max-Min)./TileSize_org);
        [xx,yy,zz]=ndgrid(0:N_tiles_new(1)-1,0:N_tiles_new(2)-1,0:N_tiles_new(3)-1);
        TilePositions_new=1+[xx(:).*TileSize_org(1),yy(:).*TileSize_org(2),zz(:).*TileSize_org(3)];
        %         TileHW=TileSize_org./2;
                        
        for t=1:prod(N_tiles_new)

            Tile_new = Tile(TilePositions_new(t,1):TilePositions_new(t,1)+TileSize_org(1)-1,...
                TilePositions_new(t,2):TilePositions_new(t,2)+TileSize_org(2)-1,...
                TilePositions_new(t,3):TilePositions_new(t,3)+TileSize_org(3)-1);
%             figure();imshow(max(im2uint8(Tile_new),[],3));

            SaveTile(usedDB,Tile_new,image_id,x,TilePositions,TilePositions_new,TileSize_org,i,t,z_level,N_tiles,SaveFolder,N_tiles_new,namecount)

            namecount = namecount + 1;
        end
    end
    
elseif strcmp(T.transform,'Rigid') || strcmp(T.transform,'Affine')
    % Transform the binding boxes of all stacks
    b=T.b';
    L=T.L;
    
    Verts=zeros(8,3,N_stacks);
    Verts_transformed=zeros(8,3,N_stacks);
    for i=1:N_stacks
        Verts(:,:,i)=ones(8,1)*StackPositions(i,:)-1+[1,1,1;StackSizes(i,1),1,1;1,StackSizes(i,2),1;StackSizes(i,1),StackSizes(i,2),1;...
            1,1,StackSizes(i,3);StackSizes(i,1),1,StackSizes(i,3);1,StackSizes(i,2),StackSizes(i,3);StackSizes(i,1),StackSizes(i,2),StackSizes(i,3)];
        Verts_transformed(:,:,i)=Verts(:,:,i)*L(:,(i-1)*3+1:i*3)'+b(i,:);
    end
    
    % Find the extent of the transformed brain space
    Max=max(max(Verts_transformed,[],3),[],1);
    Min=min(min(Verts_transformed,[],3),[],1);
    Global_shift=1-Min;
    
    % Tile this space
    N_tiles=ceil(ceil(Max-Min)./TileSize);
    [xx,yy,zz]=ndgrid(0:N_tiles(1)-1,0:N_tiles(2)-1,0:N_tiles(3)-1);
    TilePositions=1+[xx(:).*TileSize(1),yy(:).*TileSize(2),zz(:).*TileSize(3)];
    
    % Find tile vertices
    TileVerts=zeros(8,3,prod(N_tiles));
    for i=1:prod(N_tiles)
        TileVerts(:,:,i)=ones(8,1)*TilePositions(i,:)-1+[1,1,1;TileSize(1),1,1;1,TileSize(2),1;TileSize(1),TileSize(2),1;...
            1,1,TileSize(3);TileSize(1),1,TileSize(3);1,TileSize(2),TileSize(3);TileSize(1),TileSize(2),TileSize(3)];
    end
    namecount = 1;
    for i=1:prod(N_tiles)
        % Find stacks that overlap the tile after the transform
        StackInd=[];
        for j=1:N_stacks
            if Overlap(Verts_transformed(:,:,j),TileVerts(:,:,i)-ones(8,1)*Global_shift)
                StackInd=[StackInd,j];
            end
        end
        
        Tile=zeros(TileSize,StackClass);
        for j=1:length(StackInd)
            

            if params.RE.useHDF5
                X = im2uint8(hdf5read([DataFolder,'\tmp\temp_',num2str(StackInd(j)),'.h5'], '/dataset1'));
            else
                pth=StackList{StackInd(j)}(1:find(StackList{StackInd(j)}=='\',1,'last'));
                file_list=StackList{StackInd(j)}(find(StackList{StackInd(j)}=='\',1,'last')+1:end);
                X=ImportStack([pth,file_list]);
            end
            %             [X,~,~]=ImportStackJ(pth,{file_list},1,1,1,'Max');
            
            % inverse transform tile voxels
            [xx,yy,zz]=ndgrid(1:TileSize(1),1:TileSize(2),1:TileSize(3));
            xyz=[xx(:)+TilePositions(i,1)-1-Global_shift(1)-b(StackInd(j),1),yy(:)+TilePositions(i,2)-1-Global_shift(2)-b(StackInd(j),2),zz(:)+TilePositions(i,3)-1-Global_shift(3)-b(StackInd(j),3)];
            clear xx yy zz
            xyz=xyz*inv(L(:,(StackInd(j)-1)*3+1:StackInd(j)*3)');
            xyz=round(xyz-ones(prod(TileSize),1)*(StackPositions(StackInd(j),:)-1)); % coordinates within the stack
            
            ind=(xyz(:,1)>=1 & xyz(:,1)<=StackSizes(StackInd(j),1) & ...
                xyz(:,2)>=1 & xyz(:,2)<=StackSizes(StackInd(j),2) & ...
                xyz(:,3)>=1 & xyz(:,3)<=StackSizes(StackInd(j),3));
            ind2=xyz(ind,1)+(xyz(ind,2)-1).*StackSizes(StackInd(j),1)+(xyz(ind,3)-1).*(StackSizes(StackInd(j),1)*StackSizes(StackInd(j),2));
            Tile(ind)=max(Tile(ind),X(ind2));
        end
        Tile = uint8(double(Tile)./double(MaxIntensityValue)*255);
        %         Tile = uint8(double(Tile));
        
        N_tiles1 = 1;
        Verts=zeros(8,3,N_tiles1);
        
        Verts(:,:)=ones(8,1)*TilePositions(i,:)-1+[1,1,1;TileSize(1),1,1;1,TileSize(2),1;TileSize(1),TileSize(2),1;...
            1,1,TileSize(3);TileSize(1),1,TileSize(3);1,TileSize(2),TileSize(3);TileSize(1),TileSize(2),TileSize(3)];
        
        
        % Find the extent of the transformed brain space
        Max=max(max(Verts,[],3),[],1);
        Min=min(min(Verts,[],3),[],1);
        
        % Tile this space
        N_tiles_new=ceil(ceil(Max-Min)./TileSize_org);
        [xx,yy,zz]=ndgrid(0:N_tiles_new(1)-1,0:N_tiles_new(2)-1,0:N_tiles_new(3)-1);
        TilePositions_new=1+[xx(:).*TileSize_org(1),yy(:).*TileSize_org(2),zz(:).*TileSize_org(3)];
        
        
        for t=1:prod(N_tiles_new)
            %             Tile_new = zeros(TileSize_org);
            Tile_new = Tile(TilePositions_new(t,1):TilePositions_new(t,1)+TileSize_org(1)-1,...
                TilePositions_new(t,2):TilePositions_new(t,2)+TileSize_org(2)-1,...
                TilePositions_new(t,3):TilePositions_new(t,3)+TileSize_org(3)-1);
            %              figure();imshow(max(Tile_new,[],3),[0 max(Tile_new(:))]);
            
            TileName=num2str(namecount,['%0',num2str(fix(log10(prod(N_tiles)*prod(N_tiles_new)))+1),'.0f']);
            namecount = namecount + 1;
            if usedDB
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
            else
                %                 TileName=num2str(i,['%0',num2str(fix(log10(prod(N_tiles)))+1),'.0f']);
                mkdir([SaveFolder,TileName]);
                saveastiff(im2uint8(Tile), [SaveFolder,TileName,'\',TileName,'.tif']);
                disp(['Tile ',num2str(i),'/',num2str(prod(N_tiles)),' created.']);
            end
            %         [ii,jj,kk]=ind2sub(N_tiles,i);
            %         figure(kk)
            %         subplot(N_tiles(1),N_tiles(2),(ii-1)*N_tiles(2)+jj), imshow(max(Tile,[],3))
        end
        
        
        %         [ii,jj,kk]=ind2sub(N_tiles,i);
        %         figure(kk)
        %         subplot(N_tiles(1),N_tiles(2),(ii-1)*N_tiles(2)+jj), imshow(max(Tile,[],3))
        %         drawnow
    end
    
elseif strcmp(T.transform,'Nonrigid')
end
% save: TileNames TilePositions TileSizes ZoomLevels



if usedDB ==1
    close(conn)
    clear conn
    disp("Done with DB");
end


if usedDB ==1 || usedDB ==2 || usedDB ==4 
    TileInfo.N_tiles=N_tiles;
    TileInfo.TileNames=num2cell(num2str((1:prod(N_tiles))',['%0',num2str(fix(log10(prod(N_tiles)))+1),'.0f']),2);
    TileInfo.TilePositions=TilePositions;
    TileInfo.TileSizes=repmat(TileSize,prod(N_tiles),1);
    TileInfo.ZoomLevels=ones(prod(N_tiles),1);
    if ~exist(SaveFolder, 'dir')
        mkdir(SaveFolder)
    end
    save([SaveFolder,'TileInfo.mat'],'TileInfo')
    T = table(TileInfo.TileNames,TileInfo.TilePositions,TileInfo.TileSizes,TileInfo.ZoomLevels);
    writetable(T,[SaveFolder,'TileInfo.csv'],'WriteVariableNames',false)
elseif usedDB ==3
    
    json_string = ['{"data_type": "uint8", "num_channels": 1, "scales": [{"chunk_sizes": [[',num2str(TileSize_org(1)),',',num2str(TileSize_org(2)),',',num2str(TileSize_org(3)),']], "encoding": "jpeg", "key": "image", "resolution": [128, 128, 128], "size": [',num2str(TilePositions(end,2)+TileSize_org(2)),',',num2str(TilePositions(end,1)+TileSize_org(1)),',',num2str(TilePositions(end,3)+TileSize_org(3)),'], "voxel_offset": [0, 0, 0]}], "type": "image"}'];
    
%   jsonencode(containers.Map( {'data_type','num_channels','scales'}, [1,1,1]))

    fileID = fopen([SaveFolder,'/info'],'w');
    nbytes = fprintf(fileID,json_string);
    fclose(fileID);
elseif usedDB == 5
    createXML(TilePositions,TileSize_org,SaveFolder);
end



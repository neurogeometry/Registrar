% This function creates seamless image tiles from the dataset of original image stacks and stack transformations.

function Retiling(LogHandle,StackList,StackPositions,StackSizes,T,DataFolder,outputType,Seq_Par,Par_workers,DBFile)
% addpath('../Functions');
parameters
paramsREuseHDF5=paramsREuseHDF5;
paramsRERemoveBlack = paramsRERemoveBlack;
paramsBigTileSize=paramsBigTileSize;
paramsFinalTileSize=paramsFinalTileSize;
discovery = 1;

% if discovery
%     DBFile = 0;
% end
paramsEmptyVoxelsValue = paramsEmptyVoxelsValue;
options.overwrite = 1;
% InfoImage=imfinfo(char(StackList(1,1)));
% M_tiles = ceil(InfoImage(1).Height/paramsBigTileSize(1));
% paramsBigTileSize = paramsBigTileSize.*2;
Trimimage = 0;
image_id = 1;
x = [];
z_level = 1;
if Trimimage
    w = 60;
    h = 5;
    StackSizes = StackSizes - [h*2,w*2,0];
end
[filepath,~,ext] = fileparts(char(StackList(1,1)));

if paramsREuseHDF5
    X = hdf5read([DataFolder,'/tmp/temp_',num2str(1),'.h5'], '/dataset1');
    MaxIntensityValue = max(X(:));
    StackClass = class(hdf5read([DataFolder,'/tmp/temp_',num2str(1),'.h5'], '/dataset1'));
else
    if size(ext,2)>1
        InfoImage=imfinfo(char(StackList(1,1)));
        MaxIntensityValue = InfoImage(1).MaxSampleValue;
        StackClass = class(imread(char(StackList(1,1))));
    else
        allfiles = dir(filepath);
        %         InfoImage=imfinfo(char([allfiles(3).folder,'/',allfiles(3).name]));
        InfoImage=imfinfo(char([filepath,'/',allfiles(3).name]));
        if strcmp(InfoImage.Format,'jpg')
            %             FinalImage=ImportStack([allfiles(3).folder,'/'],StackSizes(1,:));
            FinalImage=ImportStack([filepath,'/'],StackSizes(1,:));
            MaxIntensityValue = max(FinalImage(:));
        else
            MaxIntensityValue = InfoImage(1).MaxSampleValue;
        end
        %         StackClass = class(imread(char([allfiles(3).folder,'/',allfiles(3).name])));
        StackClass = class(imread(char([filepath,'/',allfiles(3).name])));
    end
end


if outputType == 1
    SpecimenName = extractBefore(extractAfter(DataFolder,'MicroscopeFiles/Results-'),'_StackList');
    %     DBFile = [pwd,'/',DataFolder,'/nctracer.db'];
    %     DBFile = 'E:/TilesCreation/NCTracerWeb/New/NCtracerWeb-master/NCtracerWeb-master/NCT-Web/data/nctracer.db';
    % for connection help: https://www.mathworks.com/help/database/ug/sqlite-jdbc-windows.html#bt8kopj-1
    %     delete(DBFile);
    if exist(DBFile) > 0
        conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:',DBFile]);
        x = conn.Handle;
    else
        
        conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:',DBFile]);
        x = conn.Handle;
        sqlquery_createEDGE = 'CREATE TABLE edge (id integer primary key autoincrement not null, trace_id integer not null, v1_id integer not null, v2_id integer not null)';
        curs = exec(conn,sqlquery_createEDGE);
        close(curs);
        
        sqlquery_createIMAGE = 'CREATE TABLE image (id integer primary key autoincrement not null, name text not null)';
        curs = exec(conn,sqlquery_createIMAGE);
        close(curs);
        
        sqlquery_createPIX = 'CREATE TABLE "pix" ( `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,`image_id` integer NOT NULL, `x` integer NOT NULL, `y` integer NOT NULL, `z` integer NOT NULL, `x_dim` integer NOT NULL, `y_dim` integer NOT NULL, `z_dim` integer NOT NULL, `pixels` blob NOT NULL,`zoom_out` INTEGER NOT NULL)';
        curs = exec(conn,sqlquery_createPIX);
        close(curs);
        
        sqlquery_createTRACE = 'CREATE TABLE trace (id integer primary key autoincrement not null, name text not null, image_id integer not null)';
        curs = exec(conn,sqlquery_createTRACE);
        close(curs);
        
        sqlquery_createbatchlog = 'CREATE TABLE batchlog (id INTEGER PRIMARY KEY AUTOINCREMENT, job TEXT NOT NULL, description TEXT NOT NULL, created timestamp NOT NULL)';
        curs = exec(conn,sqlquery_createbatchlog);
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
    q2 = 'select seq from sqlite_sequence where name=''image''';
    curs = exec(conn, q2);
    curs = fetch(curs);
    image_id = curs.Data{1}
    z_level = 1;
end

SaveFolder=[DataFolder,params.RE.savefolder,'Zoom1/'];
% for test
% SaveFolder = 'E:/TilesCreation/Data/MouselightFull_Neuroglancer/';

N_stacks=length(StackList);
if outputType == 3
    mkdir([SaveFolder,'image']);
elseif outputType == 4
    mkdir([SaveFolder,'Nifti']);
elseif outputType == 5
    mkdir([SaveFolder,'HDF5']);
elseif outputType == 6
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
        %         Verts_transformed(:,:,i)=Verts(:,:,i)+b(i,:);
        Verts_transformed(:,:,i)=Verts(:,:,i)+ones(size(Verts(:,:,i),1),1)*b(i,:);
    end
    
    % Find the extent of the transformed brain space
    Max=max(max(Verts_transformed,[],3),[],1);
    Min=min(min(Verts_transformed,[],3),[],1);
    Global_shift=1-Min;
    
    StackPositionsTransformed=round(StackPositions+b+ones(N_stacks,1)*Global_shift);
    if Trimimage
        StackPositionsTransformed = StackPositionsTransformed + [h,w,0];
    end
    StackCentersTransformed=StackPositionsTransformed+(StackSizes-1)./2;
    
    % Tile this space
    N_tiles=ceil((Max-Min)./paramsBigTileSize);
    [xx,yy,zz]=ndgrid(0:N_tiles(1)-1,0:N_tiles(2)-1,0:N_tiles(3)-1);
    TilePositions=1+[xx(:).*paramsBigTileSize(1),yy(:).*paramsBigTileSize(2),zz(:).*paramsBigTileSize(3)];
    TileHW=paramsBigTileSize./2;
    
    reduction=paramsBigTileSize./paramsFinalTileSize;
    [xx,yy,zz]=ind2sub(reduction,(1:prod(reduction))');
    
    
    
    FinalTilePositions=zeros(prod(reduction)*prod(N_tiles),3);
    for i=1:prod(N_tiles)
        FinalTilePositions((i-1)*prod(reduction)+1:i*prod(reduction),:)=TilePositions(i,:)+([xx,yy,zz]-1).*(ones(prod(reduction),1)*paramsFinalTileSize);
    end
    %     remove_FinalTile_ind=(FinalTilePositions(:,1)>Max(1) | FinalTilePositions(:,2)>Max(2) | FinalTilePositions(:,3)>Max(3));
    %     FinalTile_name=(1:length(remove_FinalTile_ind))';
    %     FinalTile_name(remove_FinalTile_ind)=nan;
    %     FinalTile_name=FinalTile_name-cumsum(remove_FinalTile_ind);
    %     FinalN_tiles=N_tiles.*reduction;
    FinalTilePositions =[];
    
    if Seq_Par > 1
        if ~discovery
            parpool(Par_workers)
        end
        parfor i=1:prod(N_tiles)
            if outputType ==1
                conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:',DBFile]);
                x = conn.Handle;
            else
                x = [];
            end
            TileCenter=TilePositions(i,:)+(paramsBigTileSize-1)./2;
            StackInd=find(abs(StackCentersTransformed(:,1)-TileCenter(1))<StackHW(:,1)+TileHW(1) & abs(StackCentersTransformed(:,2)-TileCenter(2))<StackHW(:,2)+TileHW(2) & abs(StackCentersTransformed(:,3)-TileCenter(3))<StackHW(:,3)+TileHW(3));
            %            StackInd=find(abs(StackCentersTransformed(:,1)-ones(N_stacks,1)*TileCenter(1))<StackHW(:,1)+ones(N_stacks,1)*TileHW(1) & abs(StackCentersTransformed(:,2)-ones(N_stacks,1)*TileCenter(2))<StackHW(:,2)+ones(N_stacks,1)*TileHW(2) & abs(StackCentersTransformed(:,3)-ones(N_stacks,1)*TileCenter(3))<StackHW(:,3)+ones(N_stacks,1)*TileHW(3));
            Tile=zeros(paramsBigTileSize,'uint8');
            Tile_logical=false(paramsBigTileSize);
            
            
            for j=1:length(StackInd)
                if paramsREuseHDF5
                    X = hdf5read([DataFolder,'/tmp/temp_',num2str(StackInd(j)),'.h5'], '/dataset1');
                    %                     if Trimimage
                    %                         X=X(h:end-h,w:end-w,:);
                    %                     end
                else
                    %                     pth=StackList{StackInd(j)}(1:find(StackList{StackInd(j)}=='/',1,'last'));
                    %                     file_list=StackList{StackInd(j)}(find(StackList{StackInd(j)}=='/',1,'last')+1:end);
                    X=ImportStack(StackList{StackInd(j)},StackSizes(StackInd(j),:));
                    %                     if Trimimage
                    %                         X=X(h:end-h,w:end-w,:);
                    %                     end
                end
                X = uint8(double(X)./double(MaxIntensityValue)*255);
                
                TileStart=max(StackPositionsTransformed(StackInd(j),:)-TilePositions(i,:)+1,1);
                TileEnd=min(StackPositionsTransformed(StackInd(j),:)+StackSizes(StackInd(j),:)-TilePositions(i,:),paramsBigTileSize);
                StackStart=max(1,TilePositions(i,:)-StackPositionsTransformed(StackInd(j),:)+1);
                StackEnd=min(StackSizes(StackInd(j),:),TilePositions(i,:)+paramsBigTileSize-StackPositionsTransformed(StackInd(j),:));
                
                if Trimimage
                    Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=...
                        X(StackStart(1):StackEnd(1),StackStart(2):StackEnd(2),StackStart(3):StackEnd(3));
                else
                    Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=...
                        max(Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3)),...
                        X(StackStart(1):StackEnd(1),StackStart(2):StackEnd(2),StackStart(3):StackEnd(3)));
                    Tile_logical(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=true;
                end
                
            end
            
            Tile = uint8(double(Tile)./double(MaxIntensityValue)*255);
            Tile(Tile_logical==false)=paramsEmptyVoxelsValue;
            FinalTilePositions=ones(prod(reduction),1)*TilePositions(i,:)+([xx,yy,zz]-1).*(ones(prod(reduction),1)*paramsFinalTileSize);
            for jj=1:prod(reduction)
                if (FinalTilePositions(jj,1)<=Max(1) && FinalTilePositions(jj,2)<=Max(2) && FinalTilePositions(jj,3)<=Max(3))
                    [xxx,yyy,zzz]=ind2sub(reduction,jj);
                    FinalTile=Tile((xxx-1)*paramsFinalTileSize(1)+1:xxx*paramsFinalTileSize(1),(yyy-1)*paramsFinalTileSize(2)+1:yyy*paramsFinalTileSize(2),(zzz-1)*paramsFinalTileSize(3)+1:zzz*paramsFinalTileSize(3));
                    if ~isempty(find(FinalTile(:)~=paramsEmptyVoxelsValue,1,'first'))
                        SaveTile(LogHandle,outputType,FinalTile,image_id,x,FinalTilePositions(jj,:),paramsFinalTileSize,z_level,SaveFolder)
                    end
                end
            end
        end
        if ~discovery
            delete(gcp)
        end
    else
        for i=1:prod(N_tiles)
            TileCenter=TilePositions(i,:)+(paramsBigTileSize-1)./2;
            StackInd=find(abs(StackCentersTransformed(:,1)-TileCenter(1))<StackHW(:,1)+TileHW(1) & abs(StackCentersTransformed(:,2)-TileCenter(2))<StackHW(:,2)+TileHW(2) & abs(StackCentersTransformed(:,3)-TileCenter(3))<StackHW(:,3)+TileHW(3));
            %           StackInd=find(abs(StackCentersTransformed(:,1)-ones(N_stacks,1)*TileCenter(1))<StackHW(:,1)+ones(N_stacks,1)*TileHW(1) & abs(StackCentersTransformed(:,2)-ones(N_stacks,1)*TileCenter(2))<StackHW(:,2)+ones(N_stacks,1)*TileHW(2) & abs(StackCentersTransformed(:,3)-ones(N_stacks,1)*TileCenter(3))<StackHW(:,3)+ones(N_stacks,1)*TileHW(3));
            
            %             Tile=nan(paramsBigTileSize);
            %Tile=ones(paramsBigTileSize,'uint8')*paramsEmptyVoxelsValue;
            Tile=zeros(paramsBigTileSize,'uint8');
            Tile_logical=false(paramsBigTileSize);
            for j=1:length(StackInd)
                if paramsREuseHDF5
                    X = hdf5read([DataFolder,'/tmp/temp_',num2str(StackInd(j)),'.h5'], '/dataset1');
                    if Trimimage
                        X=X(h:end-h,w:end-w,:);
                    end
                else
                    %pth=StackList{StackInd(j)}(1:find(StackList{StackInd(j)}=='/',1,'last'));
                    %file_list=StackList{StackInd(j)}(find(StackList{StackInd(j)}=='/',1,'last')+1:end);
                    X=ImportStack(StackList{StackInd(j)},StackSizes(StackInd(j),:));
                    if Trimimage
                        X=X(h:end-h,w:end-w,:);
                    end
                end
                X = uint8(double(X)./double(MaxIntensityValue)*255);
                
                TileStart=max(StackPositionsTransformed(StackInd(j),:)-TilePositions(i,:)+1,1);
                TileEnd=min(StackPositionsTransformed(StackInd(j),:)+StackSizes(StackInd(j),:)-TilePositions(i,:),paramsBigTileSize);
                StackStart=max(1,TilePositions(i,:)-StackPositionsTransformed(StackInd(j),:)+1);
                StackEnd=min(StackSizes(StackInd(j),:),TilePositions(i,:)+paramsBigTileSize-StackPositionsTransformed(StackInd(j),:));
                
                if Trimimage
                    Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=...
                        X(StackStart(1):StackEnd(1),StackStart(2):StackEnd(2),StackStart(3):StackEnd(3));
                else
                    %                     temp = Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3));
                    %                     temp(temp==paramsEmptyVoxelsValue)=0;
                    %                     Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=...
                    %                         max(temp,X(StackStart(1):StackEnd(1),StackStart(2):StackEnd(2),StackStart(3):StackEnd(3)));
                    
                    Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=...
                        max(Tile(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3)),...
                        X(StackStart(1):StackEnd(1),StackStart(2):StackEnd(2),StackStart(3):StackEnd(3)));
                    Tile_logical(TileStart(1):TileEnd(1),TileStart(2):TileEnd(2),TileStart(3):TileEnd(3))=true;
                end
            end
            Tile(Tile_logical==false)=paramsEmptyVoxelsValue;
            
            FinalTilePositions=ones(prod(reduction),1)*TilePositions(i,:)+([xx,yy,zz]-1).*(ones(prod(reduction),1)*paramsFinalTileSize);
            for jj=1:prod(reduction)
                if (FinalTilePositions(jj,1)<=Max(1) && FinalTilePositions(jj,2)<=Max(2) && FinalTilePositions(jj,3)<=Max(3))
                    [xxx,yyy,zzz]=ind2sub(reduction,jj);
                    FinalTile=Tile((xxx-1)*paramsFinalTileSize(1)+1:xxx*paramsFinalTileSize(1),(yyy-1)*paramsFinalTileSize(2)+1:yyy*paramsFinalTileSize(2),(zzz-1)*paramsFinalTileSize(3)+1:zzz*paramsFinalTileSize(3));
                    if ~isempty(find(FinalTile(:)~=paramsEmptyVoxelsValue,1,'first'))
                        SaveTile(LogHandle,outputType,FinalTile,image_id,x,FinalTilePositions(jj,:),paramsFinalTileSize,z_level,SaveFolder,N_tiles)
                    end
                end
            end
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
    N_tiles=ceil(ceil(Max-Min)./paramsBigTileSize);
    [xx,yy,zz]=ndgrid(0:N_tiles(1)-1,0:N_tiles(2)-1,0:N_tiles(3)-1);
    TilePositions=1+[xx(:).*paramsBigTileSize(1),yy(:).*paramsBigTileSize(2),zz(:).*paramsBigTileSize(3)];
    
    reduction=paramsBigTileSize./paramsFinalTileSize;
    [xxxx,yyyy,zzzz]=ind2sub(reduction,(1:prod(reduction))');
    
    % Find tile vertices
    TileVerts=zeros(8,3,prod(N_tiles));
    for i=1:prod(N_tiles)
        TileVerts(:,:,i)=ones(8,1)*TilePositions(i,:)-1+[1,1,1;paramsBigTileSize(1),1,1;1,paramsBigTileSize(2),1;paramsBigTileSize(1),paramsBigTileSize(2),1;...
            1,1,paramsBigTileSize(3);paramsBigTileSize(1),1,paramsBigTileSize(3);1,paramsBigTileSize(2),paramsBigTileSize(3);paramsBigTileSize(1),paramsBigTileSize(2),paramsBigTileSize(3)];
    end
    
    [xx,yy,zz]=ndgrid(1:paramsBigTileSize(1),1:paramsBigTileSize(2),1:paramsBigTileSize(3));
    
    if Seq_Par > 1
        if ~discovery
            parpool(Par_workers)
        end
        parfor i=1:prod(N_tiles)
            % Find stacks that overlap the tile after the transform
            if outputType ==1
                conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:',DBFile]);
                x = conn.Handle;
            end
            StackInd=[];
            for j=1:N_stacks
                if Overlap(Verts_transformed(:,:,j),TileVerts(:,:,i)-ones(8,1)*Global_shift)
                    StackInd=[StackInd,j];
                end
            end
            
            if paramsRERemoveBlack
                Tile=ones(paramsBigTileSize,'uint8')*paramsEmptyVoxelsValue;
            else
                Tile=zeros(paramsBigTileSize,'uint8');
            end
            for j=1:length(StackInd)
                
                if paramsREuseHDF5
                    X = hdf5read([DataFolder,'/tmp/temp_',num2str(StackInd(j)),'.h5'], '/dataset1');
                else
                    %                     pth=StackList{StackInd(j)}(1:find(StackList{StackInd(j)}=='/',1,'last'));
                    %                     file_list=StackList{StackInd(j)}(find(StackList{StackInd(j)}=='/',1,'last')+1:end);
                    X=ImportStack(StackList{StackInd(j)},StackSizes(StackInd(j),:));
                end
                
                Linv=inv(L(:,(StackInd(j)-1)*3+1:StackInd(j)*3)');
                
                temp1=(TilePositions(i,1)-1-Global_shift(1)-b(StackInd(j),1))*Linv(1,1)+(TilePositions(i,2)-1-Global_shift(2)-b(StackInd(j),2))*Linv(2,1)+(TilePositions(i,3)-1-Global_shift(3)-b(StackInd(j),3))*Linv(3,1)-(StackPositions(StackInd(j),1)-1);
                temp2=(TilePositions(i,1)-1-Global_shift(1)-b(StackInd(j),1))*Linv(1,2)+(TilePositions(i,2)-1-Global_shift(2)-b(StackInd(j),2))*Linv(2,2)+(TilePositions(i,3)-1-Global_shift(3)-b(StackInd(j),3))*Linv(3,2)-(StackPositions(StackInd(j),2)-1);
                temp3=(TilePositions(i,1)-1-Global_shift(1)-b(StackInd(j),1))*Linv(1,3)+(TilePositions(i,2)-1-Global_shift(2)-b(StackInd(j),2))*Linv(2,3)+(TilePositions(i,3)-1-Global_shift(3)-b(StackInd(j),3))*Linv(3,3)-(StackPositions(StackInd(j),3)-1);
                
                xx2=round(xx*Linv(1,1)+yy*Linv(2,1)+zz*Linv(3,1)+temp1);
                yy2=round(xx*Linv(1,2)+yy*Linv(2,2)+zz*Linv(3,2)+temp2);
                zz2=round(xx*Linv(1,3)+yy*Linv(2,3)+zz*Linv(3,3)+temp3);
                
                ind=(xx2>=1 & xx2<=StackSizes(StackInd(j),1) & yy2>=1 & yy2<=StackSizes(StackInd(j),2) & zz2>=1 & zz2<=StackSizes(StackInd(j),3));
                ind2=xx2+(yy2-1).*StackSizes(StackInd(j),1)+(zz2-1).*(StackSizes(StackInd(j),1)*StackSizes(StackInd(j),2));
                if paramsRERemoveBlack
                    Tile(ind) = 0;
                end
                Tile(ind)=max(Tile(ind),X(ind2(ind)));
            end
            
            Tile = uint8(double(Tile)./double(MaxIntensityValue)*255);
            FinalTilePositions=TilePositions(i,:)+([xxxx,yyyy,zzzz]-1).*(ones(prod(reduction),1)*paramsFinalTileSize);
            for jj=1:prod(reduction)
                if (FinalTilePositions(jj,1)<=Max(1) && FinalTilePositions(jj,2)<=Max(2) && FinalTilePositions(jj,3)<=Max(3))
                    [xxx,yyy,zzz]=ind2sub(reduction,jj);
                    FinalTile=Tile((xxx-1)*paramsFinalTileSize(1)+1:xxx*paramsFinalTileSize(1),(yyy-1)*paramsFinalTileSize(2)+1:yyy*paramsFinalTileSize(2),(zzz-1)*paramsFinalTileSize(3)+1:zzz*paramsFinalTileSize(3));
                    if ~isempty(find(FinalTile(:),1,'first'))
                        SaveTile(LogHandle,outputType,FinalTile,image_id,x,FinalTilePositions(jj,:),paramsFinalTileSize,z_level,SaveFolder)
                    end
                end
            end
            
        end
        if ~discovery
            delete(gcp)
        end
    else
        for i=1:prod(N_tiles)
            % Find stacks that overlap the tile after the transform
            StackInd=[];
            for j=1:N_stacks
                if Overlap(Verts_transformed(:,:,j),TileVerts(:,:,i)-ones(8,1)*Global_shift)
                    StackInd=[StackInd,j];
                end
            end
            
            if paramsRERemoveBlack
                Tile=ones(paramsBigTileSize,'uint8')*paramsEmptyVoxelsValue;
            else
                Tile=zeros(paramsBigTileSize,'uint8');
            end
            for j=1:length(StackInd)
                
                if paramsREuseHDF5
                    X = hdf5read([DataFolder,'/tmp/temp_',num2str(StackInd(j)),'.h5'], '/dataset1');
                else
                    %                     pth=StackList{StackInd(j)}(1:find(StackList{StackInd(j)}=='/',1,'last'));
                    %                     file_list=StackList{StackInd(j)}(find(StackList{StackInd(j)}=='/',1,'last')+1:end);
                    X=ImportStack(StackList{StackInd(j)},StackSizes(StackInd(j),:));
                end
                
                Linv=inv(L(:,(StackInd(j)-1)*3+1:StackInd(j)*3)');
                
                temp1=(TilePositions(i,1)-1-Global_shift(1)-b(StackInd(j),1))*Linv(1,1)+(TilePositions(i,2)-1-Global_shift(2)-b(StackInd(j),2))*Linv(2,1)+(TilePositions(i,3)-1-Global_shift(3)-b(StackInd(j),3))*Linv(3,1)-(StackPositions(StackInd(j),1)-1);
                temp2=(TilePositions(i,1)-1-Global_shift(1)-b(StackInd(j),1))*Linv(1,2)+(TilePositions(i,2)-1-Global_shift(2)-b(StackInd(j),2))*Linv(2,2)+(TilePositions(i,3)-1-Global_shift(3)-b(StackInd(j),3))*Linv(3,2)-(StackPositions(StackInd(j),2)-1);
                temp3=(TilePositions(i,1)-1-Global_shift(1)-b(StackInd(j),1))*Linv(1,3)+(TilePositions(i,2)-1-Global_shift(2)-b(StackInd(j),2))*Linv(2,3)+(TilePositions(i,3)-1-Global_shift(3)-b(StackInd(j),3))*Linv(3,3)-(StackPositions(StackInd(j),3)-1);
                
                xx2=round(xx*Linv(1,1)+yy*Linv(2,1)+zz*Linv(3,1)+temp1);
                yy2=round(xx*Linv(1,2)+yy*Linv(2,2)+zz*Linv(3,2)+temp2);
                zz2=round(xx*Linv(1,3)+yy*Linv(2,3)+zz*Linv(3,3)+temp3);
                
                ind=(xx2>=1 & xx2<=StackSizes(StackInd(j),1) & yy2>=1 & yy2<=StackSizes(StackInd(j),2) & zz2>=1 & zz2<=StackSizes(StackInd(j),3));
                ind2=xx2+(yy2-1).*StackSizes(StackInd(j),1)+(zz2-1).*(StackSizes(StackInd(j),1)*StackSizes(StackInd(j),2));
                if paramsRERemoveBlack
                    Tile(ind) = 0;
                end
                %Tile(ind)=max(Tile(ind),X(ind2(ind)));
                Tile(ind)=X(ind2(ind));
            end
            
            Tile = uint8(double(Tile)./double(MaxIntensityValue)*255);
            FinalTilePositions=TilePositions(i,:)+([xxxx,yyyy,zzzz]-1).*(ones(prod(reduction),1)*paramsFinalTileSize);
            for jj=1:prod(reduction)
                if (FinalTilePositions(jj,1)<=Max(1) && FinalTilePositions(jj,2)<=Max(2) && FinalTilePositions(jj,3)<=Max(3))
                    [xxx,yyy,zzz]=ind2sub(reduction,jj);
                    FinalTile=Tile((xxx-1)*paramsFinalTileSize(1)+1:xxx*paramsFinalTileSize(1),(yyy-1)*paramsFinalTileSize(2)+1:yyy*paramsFinalTileSize(2),(zzz-1)*paramsFinalTileSize(3)+1:zzz*paramsFinalTileSize(3));
                    if ~isempty(find(FinalTile(:),1,'first'))
                        SaveTile(LogHandle,outputType,FinalTile,image_id,x,FinalTilePositions(jj,:),paramsFinalTileSize,z_level,SaveFolder)
                    end
                end
            end
        end
    end
    
elseif strcmp(T.transform,'Nonrigid')
end

if outputType ==1
    close(conn)
    clear conn
    disp('Done with DB');
end

if outputType ==1 || outputType ==2 || outputType ==4
    %     TileInfo.N_tiles=N_tiles;
    %     TileInfo.TileNames=num2cell(num2str(FinalTile_name(~remove_FinalTile_ind),['%0',num2str(fix(log10(prod(FinalN_tiles)))+1),'.0f']),2);
    %     TileInfo.TilePositions=FinalTilePositions(~remove_FinalTile_ind,:);
    %     TileInfo.paramsBigTileSizes=repmat(paramsBigTileSize,prod(N_tiles),1);
    %     TileInfo.ZoomLevels=ones(prod(N_tiles),1);
    %     if ~exist(SaveFolder, 'dir')
    %         mkdir(SaveFolder)
    %     end
    %     save([SaveFolder,'TileInfo.mat'],'TileInfo')
    %     T = table(TileInfo.TileNames,TileInfo.TilePositions,TileInfo.paramsBigTileSizes,TileInfo.ZoomLevels);
    %     writetable(T,[SaveFolder,'TileInfo.csv'],'WriteVariableNames',false)
elseif outputType ==3
    
    json_string = ['{"data_type": "uint8", "num_channels": 1, "scales": [{"chunk_sizes": [[',num2str(paramsFinalTileSize(1)),',',num2str(paramsFinalTileSize(2)),',',num2str(paramsFinalTileSize(3)),']], "encoding": "jpeg", "key": "image", "resolution": [128, 128, 128], "size": [',num2str(TilePositions(end,2)+paramsFinalTileSize(2)),',',num2str(TilePositions(end,1)+paramsFinalTileSize(1)),',',num2str(TilePositions(end,3)+paramsFinalTileSize(3)),'], "voxel_offset": [0, 0, 0]}], "type": "image"}'];
    
    %   jsonencode(containers.Map( {'data_type','num_channels','scales'}, [1,1,1]))
    
    fileID = fopen([SaveFolder,'/info'],'w');
    nbytes = fprintf(fileID,json_string);
    fclose(fileID);
elseif outputType == 5
    %     remove_FinalTile_ind=(FinalTilePositions(:,1)>Max(1) | FinalTilePositions(:,2)>Max(2) | FinalTilePositions(:,3)>Max(3));
    createXML(FinalTilePositions,paramsBigTileSize,SaveFolder);
end



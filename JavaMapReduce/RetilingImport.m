% This function creates seamless image tiles from the dataset of original image stacks and stack transformations.
function RetilingImport(StackList,StackPositions,StackSizes,DataFolder,REsavefolder,outputType,DBFile,paramsBigTileSize,paramsFinalTileSize,paramsEmptyVoxelsValue,MaxIntensityValue)
% paramsRERemoveBlack = 1;
% paramsBigTileSize = [512 512 128]; % powers of 2 only
% paramsFinalTileSize = [128*4 128*4 128]; % powers of 2 only
% paramsEmptyVoxelsValue=111;
disp('Retiling Started');
disp(size(StackList));
Trimimage = 0;
image_id = 1;
x = [];
z_level = 1;
if Trimimage
    w = 60;
    h = 5;
    StackSizes = StackSizes - [h*2,w*2,0];
end
% [filepath,~,ext] = fileparts(char(StackList(1,1)));
% MaxIntensityValue = 255;
% if size(ext,2)>1
%     InfoImage=imfinfo(char(StackList(1,1)));
%     MaxIntensityValue = InfoImage(1).MaxSampleValue;
%     StackClass = class(imread(char(StackList(1,1))));
% else
%     allfiles = dir(filepath);
%     InfoImage=imfinfo(char([filepath,'/',allfiles(3).name]));
%     if strcmp(InfoImage.Format,'jpg')
%         FinalImage=ImportStack([filepath,'/'],StackSizes(1,:));
%         MaxIntensityValue = max(FinalImage(:));
%     else
%         MaxIntensityValue = InfoImage(1).MaxSampleValue;
%     end
%     StackClass = class(imread(char([filepath,'/',allfiles(3).name])));
% end

% Use DB
if outputType == 1
    SpecimenName = extractAfter(DataFolder,'Results-');
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

SaveFolder=[DataFolder,REsavefolder,'Zoom1/'];

N_stacks = size(StackList,1);
disp(StackList);
disp(N_stacks);

StackHW=StackSizes./2;
    Verts=zeros(8,3,N_stacks);
    Verts_transformed=zeros(8,3,N_stacks);
    for i=1:N_stacks
        Verts(:,:,i)=ones(8,1)*StackPositions(i,:)-1+[1,1,1;StackSizes(i,1),1,1;1,StackSizes(i,2),1;StackSizes(i,1),StackSizes(i,2),1;...
            1,1,StackSizes(i,3);StackSizes(i,1),1,StackSizes(i,3);1,StackSizes(i,2),StackSizes(i,3);StackSizes(i,1),StackSizes(i,2),StackSizes(i,3)];
        Verts_transformed(:,:,i)=Verts(:,:,i);
    end
    
    % Find the extent of the transformed brain space
    Max=max(max(Verts_transformed,[],3),[],1);
    Min=min(min(Verts_transformed,[],3),[],1);
    Global_shift=1-Min;
    
    StackPositionsTransformed=round(StackPositions+ones(N_stacks,1)*Global_shift);
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
    
    FinalTilePositions =[];

    for i=1:prod(N_tiles)
        TileCenter=TilePositions(i,:)+(paramsBigTileSize-1)./2;
        StackInd=find(abs(StackCentersTransformed(:,1)-TileCenter(1))<StackHW(:,1)+TileHW(1) & abs(StackCentersTransformed(:,2)-TileCenter(2))<StackHW(:,2)+TileHW(2) & abs(StackCentersTransformed(:,3)-TileCenter(3))<StackHW(:,3)+TileHW(3));
        Tile=zeros(paramsBigTileSize,'uint8');
        Tile_logical=false(paramsBigTileSize);
        for j=1:length(StackInd)
            disp(StackList(StackInd(j),:));
            disp(class(StackList(StackInd(j),:)));
            disp(StackSizes(StackInd(j),:));
            disp(class(StackSizes(StackInd(j),:)));
            X=ImportStack(string(StackList(StackInd(j),:)));
%             X=ImportStack(string(StackList(StackInd(j),:)),[512,512,47]);
            if Trimimage
                X=X(h:end-h,w:end-w,:);
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
        Tile(Tile_logical==false)=paramsEmptyVoxelsValue;
        
        FinalTilePositions=ones(prod(reduction),1)*TilePositions(i,:)+([xx,yy,zz]-1).*(ones(prod(reduction),1)*paramsFinalTileSize);
        for jj=1:prod(reduction)
            if (FinalTilePositions(jj,1)<=Max(1) && FinalTilePositions(jj,2)<=Max(2) && FinalTilePositions(jj,3)<=Max(3))
                [xxx,yyy,zzz]=ind2sub(reduction,jj);
                FinalTile=Tile((xxx-1)*paramsFinalTileSize(1)+1:xxx*paramsFinalTileSize(1),(yyy-1)*paramsFinalTileSize(2)+1:yyy*paramsFinalTileSize(2),(zzz-1)*paramsFinalTileSize(3)+1:zzz*paramsFinalTileSize(3));
                if ~isempty(find(FinalTile(:)~=paramsEmptyVoxelsValue,1,'first'))
                    SaveTile(outputType,FinalTile,image_id,x,FinalTilePositions(jj,:),paramsFinalTileSize,z_level,SaveFolder,N_tiles)
                end
            end
        end
    end
if outputType ==1
    close(conn)
    clear conn
    disp('Done with DB');
end

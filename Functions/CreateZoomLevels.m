% This function creates factor of 2 zoom levels, e.g. 2,4,8,16
% Zoom is uniform in xyz. Tile sizes remain the same

function NumTiles = CreateZoomLevels(handles,LogHandle,ZoomLevel,DataFolder,OutMethod,DBFile)
addpath('../Functions');
parameters;
StackClass='uint8';
SaveFolder=[DataFolder,params.RE.savefolder];
options.overwrite = 1;
prevZoomLevel = num2str(ZoomLevel/2);
PrevZoomLevelFolder=[SaveFolder,'Zoom',prevZoomLevel,'/'];

SaveFolder=[SaveFolder,'Zoom',num2str(ZoomLevel),'/'];
if OutMethod == 1
    %     DBFile = [pwd,'/',DataFolder,'/nctracer.db'];
    %     SpecimenName = [extractBefore(extractAfter(DataFolder,"MicroscopeFiles/Results-"),'_StackList'),'_Z_',num2str(ZoomLevel)];
    
    %     DBFile = 'E:/TilesCreation/NCTracerWeb/New/NCtracerWeb-master/NCtracerWeb-master/NCT-Web/data/nctracer.db';
    conn = database('','','','org.sqlite.JDBC',['jdbc:sqlite:',DBFile]);
    conn.Message
    x = conn.Handle;
    %     setdbprefs('DataReturnFormat','cellarray');
    %     insertcommand_IMAGE = 'INSERT INTO image (name) values (?)';
    %     StatementObject = x.prepareStatement(insertcommand_IMAGE);
    %     StatementObject.setObject(1,SpecimenName);
    %     StatementObject.execute
    %     close(StatementObject);
    q2 = 'select seq from sqlite_sequence where name=''image''';
    % execute the query and retrieve the results
    curs = exec(conn, q2);
    curs = fetch(curs);
    image_id = curs.Data{1};
end



temp=log2(ZoomLevel);
if temp<=0 || temp-fix(temp)~=0
    error('ZoomLevel must equal 2, 4, 8, 16, 32, 64, 128, or 256')
end


if OutMethod == 1
    SqlStr = ['SELECT x,y,z FROM pix where image_id = ',num2str(image_id),' and zoom_out = ',prevZoomLevel];
    result = fetch(conn,SqlStr);
    % check that ZoomLevel/2 exists
    if isempty(result)
        error(['Create zoom level',num2str(ZoomLevel/2),'first'])
    end
    Tile_x=zeros(1,size(result,1));
    Tile_y=zeros(1,size(result,1));
    Tile_z=zeros(1,size(result,1));
    for i=1:size(result,1)
        Tile_x(i)=result{i,1};
        Tile_y(i)=result{i,2};
        Tile_z(i)=result{i,3};
    end
else
    % check that ZoomLevel/2 exists
    if exist(PrevZoomLevelFolder,'dir')~=7
        error(['Create zoom level',num2str(ZoomLevel/2),'first'])
    end
    Tile_names = dir(PrevZoomLevelFolder);
    Tile_x=zeros(1,length(Tile_names)-2);
    Tile_y=zeros(1,length(Tile_names)-2);
    Tile_z=zeros(1,length(Tile_names)-2);
    for i=1:length(Tile_names)-3
        temp=Tile_names(i+2).name;
        ind=find(temp=='_');
        Tile_x(i)=str2double(temp(1:ind(1)-1));
        Tile_y(i)=str2double(temp(ind(1)+1:ind(2)-1));
        Tile_z(i)=str2double(temp(ind(2)+1:end));
    end
end

N_tiles=[length(unique(Tile_x)),length(unique(Tile_y)),length(unique(Tile_z))];
N_tiles_new=ceil(N_tiles./2);

[xx,yy,zz]=ndgrid(0:N_tiles_new(1)-1,0:N_tiles_new(2)-1,0:N_tiles_new(3)-1);
NewTilePositions=1+[xx(:).*paramsFinalTileSize(1),yy(:).*paramsFinalTileSize(2),zz(:).*paramsFinalTileSize(3)];
BigTilePositions=NewTilePositions.*2-1;
NumTiles = prod(N_tiles_new);
for i=1:prod(N_tiles_new)
    Tile=paramsEmptyVoxelsValue.*ones(2*paramsFinalTileSize,StackClass);
    temp_x=BigTilePositions(i,1);
    temp_y=BigTilePositions(i,2);
    temp_z=BigTilePositions(i,3);
    temp_name=[num2str(BigTilePositions(i,1)),'_',num2str(BigTilePositions(i,2)),'_',num2str(BigTilePositions(i,3))];
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(1:paramsFinalTileSize(1),1:paramsFinalTileSize(2),1:paramsFinalTileSize(3)) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(1:paramsFinalTileSize(1),1:paramsFinalTileSize(2),1:paramsFinalTileSize(3))=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    
    temp_name=[num2str(BigTilePositions(i,1)+paramsFinalTileSize(1)),'_',num2str(BigTilePositions(i,2)),'_',num2str(BigTilePositions(i,3))];
    temp_x=BigTilePositions(i,1)+paramsFinalTileSize(1);
    temp_y=BigTilePositions(i,2);
    temp_z=BigTilePositions(i,3);
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(paramsFinalTileSize(1)+1:end,1:paramsFinalTileSize(2),1:paramsFinalTileSize(3)) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(paramsFinalTileSize(1)+1:end,1:paramsFinalTileSize(2),1:paramsFinalTileSize(3))=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    temp_x=BigTilePositions(i,1);
    temp_y=BigTilePositions(i,2)+paramsFinalTileSize(2);
    temp_z=BigTilePositions(i,3);
    temp_name=[num2str(BigTilePositions(i,1)),'_',num2str(BigTilePositions(i,2)+paramsFinalTileSize(2)),'_',num2str(BigTilePositions(i,3))];
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(1:paramsFinalTileSize(1),paramsFinalTileSize(2)+1:end,1:paramsFinalTileSize(3)) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(1:paramsFinalTileSize(1),paramsFinalTileSize(2)+1:end,1:paramsFinalTileSize(3))=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    temp_x=BigTilePositions(i,1)+paramsFinalTileSize(1);
    temp_y=BigTilePositions(i,2)+paramsFinalTileSize(2);
    temp_z=BigTilePositions(i,3);
    temp_name=[num2str(BigTilePositions(i,1)+paramsFinalTileSize(1)),'_',num2str(BigTilePositions(i,2)+paramsFinalTileSize(2)),'_',num2str(BigTilePositions(i,3))];
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(paramsFinalTileSize(1)+1:end,paramsFinalTileSize(2)+1:end,1:paramsFinalTileSize(3)) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(paramsFinalTileSize(1)+1:end,paramsFinalTileSize(2)+1:end,1:paramsFinalTileSize(3))=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    temp_x=BigTilePositions(i,1);
    temp_y=BigTilePositions(i,2);
    temp_z=BigTilePositions(i,3)+paramsFinalTileSize(3);
    temp_name=[num2str(BigTilePositions(i,1)),'_',num2str(BigTilePositions(i,2)),'_',num2str(BigTilePositions(i,3)+paramsFinalTileSize(3))];
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(1:paramsFinalTileSize(1),1:paramsFinalTileSize(2),paramsFinalTileSize(3)+1:end) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(1:paramsFinalTileSize(1),1:paramsFinalTileSize(2),paramsFinalTileSize(3)+1:end)=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    temp_x=BigTilePositions(i,1)+paramsFinalTileSize(1);
    temp_y=BigTilePositions(i,2);
    temp_z=BigTilePositions(i,3)+paramsFinalTileSize(3);
    temp_name=[num2str(BigTilePositions(i,1)+paramsFinalTileSize(1)),'_',num2str(BigTilePositions(i,2)),'_',num2str(BigTilePositions(i,3)+paramsFinalTileSize(3))];
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(paramsFinalTileSize(1)+1:end,1:paramsFinalTileSize(2),paramsFinalTileSize(3)+1:end) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(paramsFinalTileSize(1)+1:end,1:paramsFinalTileSize(2),paramsFinalTileSize(3)+1:end)=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    temp_x=BigTilePositions(i,1);
    temp_y=BigTilePositions(i,2)+paramsFinalTileSize(2);
    temp_z=BigTilePositions(i,3)+paramsFinalTileSize(3);
    temp_name=[num2str(BigTilePositions(i,1)),'_',num2str(BigTilePositions(i,2)+paramsFinalTileSize(2)),'_',num2str(BigTilePositions(i,3)+paramsFinalTileSize(3))];
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(1:paramsFinalTileSize(1),paramsFinalTileSize(2)+1:end,paramsFinalTileSize(3)+1:end) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(1:paramsFinalTileSize(1),paramsFinalTileSize(2)+1:end,paramsFinalTileSize(3)+1:end)=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    temp_x=BigTilePositions(i,1)+paramsFinalTileSize(1);
    temp_y=BigTilePositions(i,2)+paramsFinalTileSize(2);
    temp_z=BigTilePositions(i,3)+paramsFinalTileSize(3);
    temp_name=[num2str(BigTilePositions(i,1)+paramsFinalTileSize(1)),'_',num2str(BigTilePositions(i,2)+paramsFinalTileSize(2)),'_',num2str(BigTilePositions(i,3)+paramsFinalTileSize(3))];
    if OutMethod == 1
        SqlStr = ['SELECT pixels FROM pix where x = ',num2str(temp_x-1),' and y = ',num2str(temp_y-1),' and z = ',num2str(temp_z-1),' and zoom_out = ',prevZoomLevel,' and image_id = ',num2str(image_id)];
        result = fetch(conn,SqlStr);
        if ~isempty(result)
            Tile(paramsFinalTileSize(1)+1:end,paramsFinalTileSize(2)+1:end,paramsFinalTileSize(3)+1:end) = reshape(typecast(result{1},StackClass),paramsFinalTileSize(1),paramsFinalTileSize(2),paramsFinalTileSize(3));
        end
    else
        if exist([PrevZoomLevelFolder,temp_name],'dir')==7
            Tile(paramsFinalTileSize(1)+1:end,paramsFinalTileSize(2)+1:end,paramsFinalTileSize(3)+1:end)=ImportStack([PrevZoomLevelFolder,temp_name,'/',temp_name,'.tif'],paramsFinalTileSize);
        end
    end
    
    Tile=max(Tile(1:2:end-1,:,:),Tile(2:2:end,:,:));
    Tile=max(Tile(:,1:2:end-1,:),Tile(:,2:2:end,:));
    Tile=max(Tile(:,:,1:2:end-1),Tile(:,:,2:2:end));
    
    if ~isempty(find(Tile(:)~=paramsEmptyVoxelsValue,1,'first'))
        % save the tile
        TileName=[num2str(NewTilePositions(i,1)),'_',num2str(NewTilePositions(i,2)),'_',num2str(NewTilePositions(i,3))];
        if OutMethod == 1
            %             raw_im = typecast(reshape(im2uint8(Tile),1,[]),'uint8');
            raw_im = typecast(reshape(Tile,1,[]),StackClass);
            %            raw_im = reshape(Tile,1,[]);
            
            %             insertcommand = ['INSERT INTO pix (image_id,x,y,z,x_dim,y_dim,z_dim,pixels,z_level) values (?,?,?,?,?,?,?,?,?)'];
            insertcommand = 'INSERT INTO pix (image_id,x,y,z,x_dim,y_dim,z_dim,pixels,zoom_out) values (?,?,?,?,?,?,?,?,?)';
            StatementObject = x.prepareStatement(insertcommand);
            StatementObject.setObject(1,image_id);
            %             StatementObject.setObject(2,TileName);
            StatementObject.setObject(2,NewTilePositions(i,1)-1);
            StatementObject.setObject(3,NewTilePositions(i,2)-1);
            StatementObject.setObject(4,NewTilePositions(i,3)-1);
            StatementObject.setObject(5,paramsFinalTileSize(1));
            StatementObject.setObject(6,paramsFinalTileSize(2));
            StatementObject.setObject(7,paramsFinalTileSize(3));
            StatementObject.setObject(8,raw_im);
            StatementObject.setObject(9,ZoomLevel);
            StatementObject.execute;
            close(StatementObject);
            if handles.checkbox15.Value
                try
                    LogHandle.Children(2).String{end+1} = ['Tile ',TileName,' inserted.'];
                    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);drawnow
                catch
                end
            end
            disp(['Tile ',TileName,' inserted.']);
            if prod(N_tiles_new) == 1
                %                 raw_im = typecast(reshape(im2uint8(max(Tile,[],3)),1,[]),'uint8');
                raw_im = typecast(reshape(im2uint8(Tile),1,[]),'uint8');
                insertcommand = 'INSERT INTO pix (image_id,x,y,z,x_dim,y_dim,z_dim,pixels,zoom_out) values (?,?,?,?,?,?,?,?,?)';
                StatementObject = x.prepareStatement(insertcommand);
                StatementObject.setObject(1,image_id);
                %                 StatementObject.setObject(2,'Map');
                StatementObject.setObject(2,NewTilePositions(i,1)-1);
                StatementObject.setObject(3,NewTilePositions(i,2)-1);
                StatementObject.setObject(4,NewTilePositions(i,3)-1);
                StatementObject.setObject(5,paramsFinalTileSize(1));
                StatementObject.setObject(6,paramsFinalTileSize(2));
                StatementObject.setObject(7,paramsFinalTileSize(3));
                StatementObject.setObject(8,raw_im);
                StatementObject.setObject(9,ZoomLevel);
                StatementObject.execute;
                close(StatementObject);
            end
        else
            TileName=[num2str(NewTilePositions(i,1)),'_',num2str(NewTilePositions(i,2)),'_',num2str(NewTilePositions(i,3))];
            mkdir([SaveFolder,TileName]);
            saveastiff(Tile, [SaveFolder,TileName,'/',TileName,'.tif'],options);
            
            % for save as JPEG (separated z) for Joe
            %             for j=1:size(Tile,3)
            %                 TileName1=[num2str(NewTilePositions(i,1)),'_',num2str(NewTilePositions(i,1)),'_',num2str(NewTilePositions(i,1)+j-1)];
            %                 imwrite(Tile(:,:,j),[SaveFolder,TileName,'/',TileName1,'.jpg'],'jpg');
            %             end
            TileName=[num2str(NewTilePositions(i,2)),'_',num2str(NewTilePositions(i,1)),'_',num2str(NewTilePositions(i,3))];
            TileName1 = ['PNG/',TileName];
            mkdir([SaveFolder,TileName1]);
            for z = 0:size(Tile,3)-1       
                TileName_1 = [sprintf('%04d', z),'.png'];
                imwrite(Tile(:,:,z+1),[SaveFolder,TileName1,'/',TileName_1],'png');
            end
            
%             % for save as JPEG (Like Neuroglancer) for Joe
%             C1 = permute(Tile,[1 3 2]);
%             Tile_glancer1 = reshape(C1,[],size(Tile,2),1);
%             xstart1 = NewTilePositions(i,2)-1;
%             xend1 = xstart1+size(Tile,1);
%             ystart1 = NewTilePositions(i,1)-1;
%             yend1 = ystart1+size(Tile,2);
%             zstart1 = NewTilePositions(i,3)-1;
%             zend1 = zstart1+size(Tile,3);
%             TileName1 = [num2str(xstart1),'-',num2str(xend1),'_',num2str(ystart1),'-',num2str(yend1),'_',num2str(zstart1),'-',num2str(zend1)];
%             %                 saveastiff(im2uint8(Tile_glancer), [SaveFolder,'image/',TileName]);
%             %             SaveFolder1 = ['C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-Neocortical2_StackList\forJoe\Zoom',num2str(ZoomLevel),'\'];
%             imwrite(Tile_glancer1,[SaveFolder,TileName1,'.jpg'],'jpg');
            
            if handles.checkbox15.Value
                try
                    LogHandle.Children(2).String{end+1} = ['Tile ',TileName,' created.'];
                    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);drawnow
                catch
                end
            end
            disp(['Tile ',TileName,' created.']);
        end
    end

end
if OutMethod == 1
    close(conn)
    clear conn
end

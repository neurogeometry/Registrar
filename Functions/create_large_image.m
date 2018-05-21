function Tile3D = create_large_image(StackPositions_Registered,StackSizes_pixels,StackList,DataFolder,stackID)

addpath('../Functions');


save_path = DataFolder;%'../data/RegisteredTiles/';

if ispc
    disp('Running on Windows!');
    delete ../data/RegisteredTiles/*
else
    disp('Running on Linix!');
    system('rm -rf RegisteredTiles/*');
end

Neighbors = [stackID; findStackNeighbors(stackID,StackSizes_pixels,StackPositions_Registered)];
PlotTiles = 0;


tileSize = [3100 3100 260]; % in pixels

StackPositions_pixels(:,1) = max(StackPositions_Registered(:,1))-StackPositions_Registered(:,1);
StackPositions_pixels(:,2) = max(StackPositions_Registered(:,2))-StackPositions_Registered(:,2);
StackPositions_pixels(:,3) = max(StackPositions_Registered(:,2))-StackPositions_Registered(:,3);

MinBox_pixels=min(StackPositions_pixels,[],1); % box origin in pixels
MaxBox_pixels=max(StackPositions_pixels,[],1)+StackSizes_pixels(1,1:3); % box end in pixels

N_tiles=ceil((MaxBox_pixels-MinBox_pixels+1)./tileSize); % tile numbers in x, y,and z
[ii,jj,kk]=ind2sub(N_tiles,(1:prod(N_tiles))');
TilePositions_pixels=([ii,jj,kk]-1).*tileSize+MinBox_pixels; % global tile positions in pixels
%save('../data/TilePositions_pixels.mat','TilePositions_pixels');
%save('../data/StackPositions_pixels.mat','StackPositions_pixels');
% regcsv = zeros(size(TilePositions_pixels,1),7);
n_file_digits=fix(log10(prod(N_tiles)))+1;

tic
for i = size(Neighbors,1):-1:1
    try
        SourceID = Neighbors(i);
%         tifFile_1 = dir([char(StackList(SourceID,2)) '/*.tif' ]);
%         Source_Stack_Folder = [tifFile_1(1).folder,'/'];
%         Source_Stack_File = tifFile_1(1).name;
         disp(['Reading ',char(StackList(SourceID,1))]);

        IM_Source = ImportStack(char(StackList(SourceID,1)));

        if PlotTiles == 1
            DrawCube([StackPositions_pixels(SourceID,1),StackPositions_pixels(SourceID,2),StackPositions_pixels(SourceID,3)], StackSizes_pixels(SourceID,1:3));
            drawnow();
        end
        
        CS=StackPositions_pixels(SourceID,:)+(StackSizes_pixels(SourceID,1:3)-1)./2; %stack center
        CT=TilePositions_pixels+(tileSize-1)./2; %tile center
        inStackTilesIdx=find(sum(abs(CS-CT)<(StackSizes_pixels(SourceID,1:3)+tileSize)./2,2)==3)
           
        for tile =1 :length(inStackTilesIdx)%length(inStackTilesIdx):-1:1
            TileNumber = inStackTilesIdx(tile);
            Tile_File =[num2str(TileNumber,['%0.',num2str(n_file_digits),'d']),'.tif'];

            
            if exist([save_path,'\',Tile_File], 'file') == 2
                disp(['Reading Tile',save_path,'\',Tile_File]);
                Tile3D=ImportStack([save_path,'\',Tile_File]);
                options.overwrite = true;
            else
                Tile3D=zeros(tileSize,'uint16');
                options.overwrite = false;
            end
            
            if PlotTiles == 1
                figure(1);
                hold on;DrawCube([TilePositions_pixels(inStackTilesIdx(tile),1),TilePositions_pixels(inStackTilesIdx(tile),2),TilePositions_pixels(inStackTilesIdx(tile),3)], tileSize);
                drawnow();
            end
            
            TileStart = StackPositions_pixels(SourceID,:)-TilePositions_pixels(TileNumber,:) +1;
            TileStart(TileStart<0)=1;
            StackStart = TilePositions_pixels(TileNumber,:)-StackPositions_pixels(SourceID,:)+1;
            StackStart(StackStart<0)=1;
            TileEnd = (StackPositions_pixels(SourceID,:)+StackSizes_pixels(SourceID,1:3)) - (TilePositions_pixels(TileNumber,:)+tileSize);%oposite

            if TileEnd(1) >= 0 TileEnd(1)=tileSize(1); end;
            if TileEnd(2) >= 0 TileEnd(2)=tileSize(2); end;
            if TileEnd(3) >= 0 TileEnd(3)=tileSize(3); end;
            
            if TileEnd(1) <= 0 TileEnd(1)=tileSize(1)+TileEnd(1); end;
            if TileEnd(2) <= 0 TileEnd(2)=tileSize(2)+TileEnd(2); end;
            if TileEnd(3) <= 0 TileEnd(3)=tileSize(3)+TileEnd(3); end;
            
            StackEnd = StackStart + (TileEnd - TileStart);
            
            Tile3D(TileStart(2):TileEnd(2),TileStart(1):TileEnd(1),TileStart(3):TileEnd(3))= ...
                max(Tile3D(TileStart(2):TileEnd(2),TileStart(1):TileEnd(1),TileStart(3):TileEnd(3)),...
                IM_Source(StackStart(2):StackEnd(2),StackStart(1):StackEnd(1),StackStart(3):StackEnd(3)));
            options.message = false;
              saveastiff(Tile3D, [save_path,'\',Tile_File],options);
%               figure(2);imshow(max(Tile3D,[],3),[0 max(Tile3D(:))]);
%               figure(3);imshow(max(IM_Source(StackStart(2):StackEnd(2),StackStart(1):StackEnd(1),StackStart(3):StackEnd(3)),[],3),[0 max(Tile3D(:))]);
    
        end   
    catch
       disp(['Error Reading the file ',SourceID]);
    end
    
end
end
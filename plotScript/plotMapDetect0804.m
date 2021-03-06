% This scrip load detection and to make mapping
clear all
close all
path2Detect = '/Users/Administrator/Documents/MPO/BRing/Data/detections/';

ID{1} = 'MLB';
ID{2} = 'PRC';

ii=1;
file{ii} = ['DETECT_NARW_' ID{ii} '_PR_V1.mat'];
D(ii) = load([path2Detect file{ii}]);

ii=2
file{ii} = ['DETECT_NARW_' ID{ii} '_PR_V1.mat'];
D(2) = load([path2Detect file{2}]);

% Dafautl parameter
invInd = [2 1];
tDelayMax =  15;
% Modify structure
%save([path2Detect file{ii}],'-struct','D','dtime','angleAD','angleMD','energiedB','file')

%% Get the file name
for ii=1:2
    % Get file name
    tmpF = getWavName(D(ii).dtime, getDirectory(['toaster' ID{ii} ]));
    D(ii).file = tmpF;
    D(ii).nbF = numel(D(ii).file);
end
%% Find the matching time

% Index of structure to compare
indS = [1 2; 2 1];
for ii=1:2
    for jj=1:numel(D(indS(ii,1)).dtime)
        % Compare file name
        diffT = abs((D(indS(ii,2)).dtime - D(indS(ii,1)).dtime(jj)));
        diffTBol = diffT <= seconds(tDelayMax);
        if  any(diffTBol) % We have close time
            if sum(diffT) == 1   % We ony have one close time
                ind = find(diffTBol, 1,'first');
                D(indS(ii,1)).match(jj) =  ind;
            else
                [~, ind] =  min(diffT);
                D(indS(ii,1)).match(jj) =  ind;
            end
        else
            D(indS(ii,1)).match(jj) = nan;
        end
    end
end
%% Macro adjustment
%D(1).angleMD(23) =  111;
%D(2).angleMD( D(1).match(23) ) = 110



%% Creating position of all detection
ii=1;

for kk=2:-1:1
    kk;
    bolM= ~isnan(D(kk).match);
    indM = find(bolM);
    nbM = length(indM);
    
    locD = nan(nbM,2);
    distD = nan(nbM,2);
    
    for ii=1:nbM
        %ii
        if kk==1
            [tmpLoc tmpDist]= beamCrossing(ID{1},ID{2}, D(1).angleMD(indM(ii)), D(2).angleMD( D(1).match(indM(ii)) ) ,'showFig', 0);
        else
            [tmpLoc tmpDist]= beamCrossing(ID{1},ID{2}, D(1).angleMD( D(2).match(indM(ii)) ), D(2).angleMD( indM(ii) ) ,'showFig', 0);
        end
        
        if isempty(tmpLoc)
            locD(ii,:) = [nan nan];
            distD(ii,:) = [nan nan];
        else
            locD(ii,:) = tmpLoc;
            distD(ii,:) = tmpDist;
        end
    end
    D(kk).lat = nan(D(kk).nbF,1);
    D(kk).lon = nan(D(kk).nbF,1);
    D(kk).dist = nan(D(kk).nbF,1);
    
    D(kk).bolDoubleD = bolM';
    D(kk).lat(bolM) = locD(:,1);
    D(kk).lon(bolM) = locD(:,2);
    D(kk).dist(bolM) = mean(distD,2);
    
    
end

%% Saving location to .xcel file

%%
for ii=1:2
        file2 = [];
        time2 =[];
        bol2 = [];
        azimut2=[];
    for kk=1:length(D(ii).file)

        if isnan(D(ii).match(kk))
            bol2(kk) = 0;
            file2{kk} = '';
            time2{kk} = '';
            azimut2(kk) = NaN;
        else
            bol2(kk) = 1;
            file2{kk} = D(invInd(ii)).file{D(ii).match(kk)} ;
            time2{kk} = D(invInd(ii)).dtime(D(ii).match(kk)) ;
            azimut2(kk) = D(invInd(ii)).angleMD(D(ii).match(kk)) ;
        end
        
    end
    
    D(ii).file2 = file2';
    D(ii).time2 = time2';
    D(ii).azimut2 = azimut2';
end
%%
ii=1
tab1 = table(D(ii).file,D(ii).dtime',D(ii).energiedB,D(ii).angleMD,D(ii).bolDoubleD,D(ii).dist,D(ii).lat,D(ii).lon,D(ii).file2,D(ii).time2,D(ii).azimut2)
ii=2
tab2 = table(D(ii).file,D(ii).dtime',D(ii).energiedB,D(ii).angleMD,D(ii).bolDoubleD,D(ii).dist,D(ii).lat,D(ii).lon,D(ii).file2,D(ii).time2,D(ii).azimut2)
%%
tab1.Properties.VariableNames = {'file','time','energie dB','azimut', 'Find2', 'Dist', 'Lat','Lon','File2','time2','azimut2'};
tab2.Properties.VariableNames = {'file','time','energie dB','azimut', 'Find2', 'Dist', 'Lat','Lon','File2','time2','azimut2'};
%%
writetable(tab1,[path2Detect 'mlbDetectResume.xlsx'])
writetable(tab2,[path2Detect 'prcDetectResume.xlsx'])


%% Mapping detection

pos1 = getArrInfo(ID{1});
pos2 = getArrInfo(ID{2});

% Zone de la carte
lon_min = -64.4;
lon_max = -61;
lat_min =  47.5;
lat_max =  48.7;

close all
disp('Start mapping')
% Load map
gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
gebco = load([getDirectory('map') 'gebco_colormap.dat']);

% figure info
fig = figure('Units', 'normalized', 'Position', [-1,0,0.8,1]);
set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen

% Les isobathes ? contourer
isobath = [0:10:150];

% Start mapping
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');
hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
m_gshhs_f('patch',[.7 .7 .7]); % coastlines

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco);

% Information graphique
c = colorbar('location','eastoutside');
cpos = c.Position;
set(c, 'FontSize', 12, 'FontName','times')
ylabel(c,'Depth(m)')

% Hydrophone location
% ------------ Plot points ----------------
hc1 = m_plot(pos1(2),pos1(1),'ko','MarkerFaceColor','k','MarkerSize',7);
hc2 = m_plot(pos2(2),pos2(1),'ko','MarkerFaceColor','k','MarkerSize',7);

% Plot detection
for ii=1:nbM
    if ~isnan(locD(ii,1))
        m_plot(locD(ii,2),locD(ii,1),'ro','MarkerFaceColor','r','MarkerSize',7)
    end
end

print([getDirectory('fout') '/generalFigures/' 'mappingDectections_210804_wide.png'],'-dpng','-r150')
%% Plot more zoom

pos1 = getArrInfo(ID{1});
pos2 = getArrInfo(ID{2});

% Zone de la carte
lon_min = -64.4;
lon_max = -64;
lat_min =  48.5;
lat_max =  48.7;

close all
disp('Start mapping')
% Load map
gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
gebco = load([getDirectory('map') 'gebco_colormap.dat']);

% figure info
fig = figure('Units', 'normalized', 'Position', [-1,0,0.8,1]);
set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen

% Les isobathes ? contourer
isobath = [0:10:150];

% Start mapping
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');
hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
m_gshhs_f('patch',[.7 .7 .7]); % coastlines

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco);

% Information graphique
c = colorbar('location','eastoutside');
cpos = c.Position;
set(c, 'FontSize', 12, 'FontName','times')
ylabel(c,'Depth(m)')

% Hydrophone location
% ------------ Plot points ----------------
hc1 = m_plot(pos1(2),pos1(1),'ko','MarkerFaceColor','k','MarkerSize',7);
hc2 = m_plot(pos2(2),pos2(1),'ko','MarkerFaceColor','k','MarkerSize',7);

%%
% Plot detection
for ii=1:nbM
    if ~isnan(locD(ii,1))
        m_plot(locD(ii,2),locD(ii,1),'ro','MarkerFaceColor','r','MarkerSize',7)
    end
end

print([getDirectory('fout') '/generalFigures/' 'mappingDectections_210804_close.png'],'-dpng','-r150')
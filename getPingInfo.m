function [time, ploc] = getPingInfo(dataID,varargin)
% [ploc, ptime ]  = getPingLoc('cav')
% This fonction return the location and time of any set of ping emission
% For circular ping, the ping number 1 is located the most north
% Times are in UTC.
% dataID option :  aav = cricle around Anse-a-Valleau
%                  cld = Cicle around Cloridorme
%                  mlb =   "  " Malbay / prc = " " perce
%                  perpH = perpendicular transect for Honguedo 
%                  more to come....
%  Orientation : 
%   AAV = Horaire and Not reverse
%   CLD = Horaire and reverse
%   PRC = Horaire and reverse
%   MLB = Horaire and not sure the side
%   
%   last update 14/09/2021

showFig = [];
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'showfig'
                showFig = varargin{2}; 
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end
%%
switch lower(dataID)
    case 'aav'   % -------- Circle at Anse-a-Valleau -----------
        ploc = [ 49.0953 64.5370;
            49.0945 64.5328;
            49.0922 64.5303;
            49.0893 64.5304;
            49.0871 64.5329;
            49.0863 64.5370;
            49.0872 64.5410;
            49.0894 64.5435;
            49.0922 64.5434;
            49.0945 64.5409];
        time=[ datetime(2021,07,15,18,10,58);
            datetime(2021,07,15,18,13,38);
            datetime(2021,07,15,18,16,10); % 16:20
            datetime(2021,07,15,18,19,04);
            datetime(2021,07,15,18,21,56);  % Was written 18:31:37 in logbook
            datetime(2021,07,15,18,25,23);
            datetime(2021,07,15,18,28,04);
            datetime(2021,07,15,18,30,25);
            datetime(2021,07,15,18,06,03);
            datetime(2021,07,15,18,08,27);];
    case 'cld'   % -------- Circle at Cloridorme -----------
        ploc = [49.1968 64.8351;
            49.1977 64.8315;
            49.1967 64.8274;
            49.1945 64.8248;
            49.1916 64.8247;
            49.1895 64.8273;
            49.1887 64.8312;
            49.1895 64.8355;
            49.1917 64.8381;
            49.1945 64.8378;];
            % 49.1931 64.8314                   % Center array emission
        time=[ datetime(2021,07,15,14,25,05);      % Two emission at same locaion
            datetime(2021,07,15,14,29,08);
            datetime(2021,07,15,14,31,55);
            datetime(2021,07,15,14,34,07);
            datetime(2021,07,15,14,36,42);  % Also a 39:31
            datetime(2021,07,15,14,39,31); 
            datetime(2021,07,15,14,42,44);  % Alos 42:44
            datetime(2021,07,15,14,46,10);
            datetime(2021,07,15,14,19,10);
            datetime(2021,07,15,14,22,41);
            %datetime(2021,07,15,14,51,20);     % Center array emission         
            ];
        case 'mlb'   % -------- Circle at Malbay -----------
        % Malbay ping have been made by an knock on the boat
        ploc = [48.6084 64.1863;
                48.6076 64.1824;
                48.6053 64.1800;
                48.6025 64.1966;
                48.6002 64.1823;
                48.5994 64.1863;
                48.6003 64.1903;
                48.6026 64.1928;
                48.6054 64.1926;
                48.6076 64.1903;]; % center emission
                %48.6039 64.1862;];
        time=[ datetime(2021,08,01,14,30,08);
            datetime(2021,08,01,14,32,05);
            datetime(2021,08,01,14,34,18);
            datetime(2021,08,01,14,36,23);
            datetime(2021,08,01,14,38,11);
            datetime(2021,08,01,14,40,01);
            datetime(2021,08,01,14,41,43);
            datetime(2021,08,01,14,43,37);
            datetime(2021,08,01,14,45,05);
            datetime(2021,08,01,14,48,59);];
            %datetime(2021,08,01,14,49,05);  % center emission
        case 'prc'   % -------- Circle at Perce -----------
        ploc = [48.5356 64.1984;
                48.5348 64.1945;
                48.5326 64.1921;
                48.5298 64.1919;
                48.5275 64.1945;
                48.5267 64.1986;
                48.5277 64.2024;
                48.5298 64.2049;
                48.5325 64.2048;
                48.5349 64.2024;];
                %48.5473 64.1991];
        time=[ datetime(2021,08,01,13,40,37);      % Two emission at same locaion
               datetime(2021,08,01,13,43,15);
               datetime(2021,08,01,13,45,54);
               datetime(2021,08,01,13,47,38);
               datetime(2021,08,01,13,49,34);
               datetime(2021,08,01,13,51,03);
               datetime(2021,08,01,13,53,06);
               datetime(2021,08,01,13,53,01);
               datetime(2021,08,01,13,56,15);
               datetime(2021,08,01,13,58,01);];
               %datetime(2021,08,01,14,00,45);
               
               case 'perpH'   % -------- perp Honguedo -----------
        ploc = [ 49.0953 64.5370;
            49.0945 64.5328;
            49.0922 64.5303;
            49.0893 64.5304;
            49.0871 64.5329;
            49.0863 64.5370;
            49.0872 64.5410;
            49.0894 64.5435;
            49.0922 64.5434;
            49.0945 64.5409];
        time=[ datetime(2021,07,15,18,11,00);
            datetime(2021,07,15,18,14,00);
            datetime(2021,07,15,18,16,20);
            datetime(2021,07,15,18,18,58);
            datetime(2021,07,15,18,21,37);  % Was written 18:31:37 in logbook
            datetime(2021,07,15,18,25,15);
            datetime(2021,07,15,18,27,56);
            datetime(2021,07,15,18,30,12);
            datetime(2021,07,15,18,05,02);
            datetime(2021,07,15,18,08,44);];
    otherwise
        error(['Can''t location ID. Refer to heater.' ])
end

if any(showFig ==1)
    
    
        
disp('Printing map')
    try
        % Read bathymetry
        load Documents/MATLAB/Oceanographie/m_map/mmap_ex/GSL_bathy_500m.mat;
        load Documents/MATLAB/Oceanographie/m_map/mmap_ex/gebco_colormap.dat;
    catch
        try
            % Read bathymetry
            load GSL_bathy_500m.mat;
            load gebco_colormap.dat;
        catch
            error("Cant not load bathy file.")
        end
    end
    
    
% Default parameter
winSize = 0.01;

lon = -ploc(:,2);
lat = ploc(:,1);

% Zone de la carte
lon_min = min(lon)-winSize;
lon_max = max(lon)+winSize;
lat_min =  min(lat)-winSize;
lat_max =  max(lat)+winSize;

% Les isobathes à contourer
isobath = [0:10:450];

% Dimensions de la figure sur papier (en cm) 
paperwidth  = 12;
paperheight = 12;

figure(1)
clf

m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('box','fancy','fontsize',8,'linestyle','none','fontname','times');

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco_colormap);

hold on
m_contourf(LON,LAT,Z,isobath,'linestyle','none'); 
% Information graphique
c = colorbar('location','eastoutside');
%set(c,'Position',[0.2 0.2 0.4 0.02],'FontSize', 8, 'FontName','times','xAxisLocation','top');
cpos = c.Position;
set(c, 'FontSize', 8, 'FontName','times')
ylabel(c,'Depth(m)')
%m_text(-69.7,46.9,'Depth (m)','FontSize', 8, 'FontName','times','HorizontalAlignment','center');
xlabel('Longitude', 'FontSize', 10, 'FontName','times')
ylabel('Latitude', 'FontSize', 10, 'FontName','times')
set (gca,'TicklabelInterpreter','None')

m_gshhs_f('patch',[.7 .7 .7]); % coastlines
  

% ------------ Plot position ----------------

%ha=   m_plot(arrLoc(2),arrLoc(1),'*','MarkerEdgeColor','k',...
%    'MarkerFaceColor','k','MarkerSize',6) ;

hold on

% Boat location
hb=   m_plot(lon(:),lat(:),'o','MarkerEdgeColor','b',...
    'MarkerFaceColor','b','MarkerSize',4) ;

for i=1:length(lon)
m_text(lon(i)+ 0.0001 ,lat(i,1)+ 0.0001,num2str(i),'color','b');
end

%print('-dpng', '-r300','figures/WBR_map_honguedo_transNCirc.png');    
disp('Done!')
    
end


function [mLim ] = getMapLim(arrID, dist, varargin)
% Set the limit for map figure 
% To upgrade. 

mode = 'fixe';

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'mode'
            mode = varargin{2};
        otherwise
            error(['Can''t understand property: ' varargin{1}])
    end
    varargin(1:2)=[];
end



if size(dist) ==2    
    posArr = getArrInfo(arrID1);
    [dist,~,~] = m_idist(posArr(2),posArr(1),dist(2),dist(1));
end

if strcmp(mode, 'fixe')
if strcmp(arrID,'PRC') || strcmp(arrID,'MLB')
    ok=1
    if dist < 20 * 1000
        lon_min = -64.4;
        lon_max = -63.9;
        lat_min =  48.4;
        lat_max =  48.7;
    elseif dist < 80 * 1000
        lon_min = -64.4;
        lon_max = -63;
        lat_min =  48.2;
        lat_max =  48.7;
    else
        lon_min = -64.4;
        lon_max = -61;
        lat_min =  47.5;
        lat_max =  49.5;
    end
elseif strcmp(arrID,'AAV') || strcmp(arrID,'CLD')
    ok=2
    lon_min = -67;
    lon_max = -64;
    lat_min =  48.8;
    lat_max =  49.5;
else
    error('Wrong arrID.')
end

else
    posArr = getArrInfo(arrID1);
    lon_min = -max(abs([posArr(:,2) dist(:,2)])) - winSize;
    lon_max = -min(abs([posArr(:,2) dist(:,2)])) + winSize ;
    lat_min =  min([posArr(:,1) dist(:,1)]) -  winSize/2;
    lat_max =  max([posArr(:,1) dis(:,1)]) + winSize/2;
end

% Send valus to a structure to 
mLim.lonMin = lon_min; 
mLim.lonMax = lon_max;
mLim.latMin = lat_min;
mLim.latMax = lat_max;
end

function hLoc = createHLoc(cLoc,varargin)
% This script produce the théoreotical location(lat/lon) circular array
% using a center point

% Defualt parameter
theta = 0;  % Off-set from true north
nbH = 10;   % Number of hydrohpone
r = 10;     % radius of the array

%% Input argument
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'theta'
                theta = varargin{2};
            otherwise
                error('Can''t understand property in convLoc.m')
        end
        varargin(1:2)=[];
end

%% Create new loc
hLoc = nan(nbH,2);

for ii=1:nbH
    [lon,lat, phi] = m_fdist(cLoc(2),cLoc(1), ((ii-1)* 360 / nbH ) + theta , r);
    hLoc(ii,1) = lat;
    hLoc(ii,2) = lon;
end


end

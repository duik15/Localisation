function outLoc = convLoc(inLoc, varargin)
% Convert lat/lon location in a newformat

% Default value
inFrmt = 'degminsecdec';
outFrmt  ='degdec';

% Add varargin parameter here
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'infrmt'
                inFrmt = varargin{2};
            otherwise
                error('Can''t understand property in convLoc.m')
        end
        varargin(1:2)=[];
end

if strcmp(inFrmt, 'degminsecdec')
    outLoc = inLoc(1) + inLoc(2)/60 + inLoc(3) /3600;
elseif strcmp(inFrmt, 'degmindec') 
    outLoc = inLoc(1) + inLoc(2)/60 ;
else
    error('Do not understand inpout format in convLoc.m')
end

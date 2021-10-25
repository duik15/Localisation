function dirPath = getDirectory(id,varargin)


machineID = getMachine;
%%
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'machineID'
                machineID = varargin{2}; 
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end

switch lower(id)
    case 'boat'
        switch lower(machineID)
            case 'dellkev'
                dirPath = 'C:\Users\duquettek\Documents\BRing\Data\boatTrack' ;
            case 'mac'
                dirPath = '/Users/Administrator/Documents/MPO/BRing/Data/boatTrack/';s
        end
    otherwise 
        error('Can''t understand the ID in getDirectory. See help getDirectory.')
end
end
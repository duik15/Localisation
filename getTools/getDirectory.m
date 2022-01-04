function dirPath = getDirectory(id,varargin)


machineID = getMachine;
%%
while ~isempty(varargin)
    machineID = varargin{1};
    %switch lower(varargin{1})
    %    case 'machineID'
    %        machineID = varargin{2};
    %    otherwise
    %        error(['Can''t understand property: ' varargin{1}])
    %end
    %varargin(1:2)=[];
end
%%

switch lower(id)
    case 'fout'
        switch lower(machineID)
            case 'labo'
                dirPath = 'Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\';
            case 'dellkev'
                dirPath = 'C:\Users\duquettek\Documents\BRing\results\' ;
            case 'mac'
                dirPath = '/Users/Administrator/Documents/MPO/BRing/Data/results/';
        end
    case 'boat'
        switch lower(machineID)
            case 'dellkev'
                dirPath = 'C:\Users\duquettek\Documents\BRing\Data\boatTrack' ;
            case 'mac'
                dirPath = '/Users/Administrator/Documents/MPO/BRing/Data/boatTrack/';
        end
    case 'aav'
        switch lower(machineID)
            case 'dellkev'
                dirPath = 'E:\Bring_Dep_1\' ;
            case 'mac'
                dirPath = '/Users/Administrator/Documents/MPO/BRing/Data/wav/';
            case 'labo'
                dirPath = ['\\169.254.116.24\usbshare2-2\Bring_Dep_1\'];
        end
    case 'cld'
        switch lower(machineID)
            case 'dellkev'
                dirPath = 'F:\Bring_Dep_1\' ;
            case 'mac'
                dirPath = '/Users/Administrator/Documents/MPO/BRing/Data/wav/';
            case 'labo'
                dirPath = ['\\169.254.47.215\usbshare1-2\Bring_Dep_1\'];
           case 'toaster'
                dirPath = '/Volumes/BringDD4/Bring_Dep_1/CLD/';
        end
    case 'mlb'
        switch lower(machineID)
            case 'dellkev'
                dirPath = 'E:\Bring_Dep_2\' ;
            case 'mac'
                dirPath = '/Users/Administrator/Documents/MPO/BRing/Data/wav/';
            case 'labo'
                dirPath = ['\\169.254.116.24\usbshare2-2\Bring_Dep_2\'];
            case 'toaster'
                dirPath = '/Volumes/BringDD3/Bring_Dep_2/MLB/';
        end
    case 'prc'
        switch lower(machineID)
            case 'dellkev'
                dirPath = 'F:\Bring_Dep_2\' ;
            case 'mac'
                dirPath = '/Users/Administrator/Documents/MPO/BRing/Data/wav/';
            case 'labo'
                dirPath = ['\\169.254.47.215\usbshare1-2\Bring_Dep_2\'];                
        end
    case 'toasterprc'
        dirPath = '/Volumes/BringDD4/Bring_Dep_2/PRC/';
    case 'toastermlb'
        dirPath = '/Volumes/BringDD3/Bring_Dep_2/MLB/';
    case 'toastercld'
        dirPath = '/Volumes/BringDD4/Bring_Dep_1/CLD/';
    case 'toasteraav'
        dirPath = '/Volumes/BringDD3/Bring_Dep_1/AAV/';
    case 'map'
        switch lower(machineID)
            case 'dellkev'
                dirPath = 'T:\ALGORITHMES\mapTools\' ;
            case 'mac'
                dirPath = 'Documents/MATLAB/Oceanographie/m_map/mmap_ex/';
            case 'labo'
                %dirPath = ['\\169.254.47.215\usbshare1-2\Bring_Dep_2\']
        end
otherwise
    error('Can''t understand the ID in getDirectory. See help getDirectory.')
end
end
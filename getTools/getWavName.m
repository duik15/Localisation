function [fileName wi] = getWavName(dateIn, folderIn)
%This load your wave file by specifying the folder, date and time
%fileName = getWavName(datetime(2021,07,15,14,37,00),mypath)
%
%dateIn

%Loading folderIn and files informations
dirInfo = dir(folderIn);
dirInfo([dirInfo.isdir]) = [];
fileList = {dirInfo.name};

typeHL = 'LF';

% keep only the wav file
i2erase = [];
for i=1:numel(fileList)
    splitName = split(fileList{i}, '_');
    if strcmp(fileList{i}(end-2:end),'wav') == 0
        i2erase = [i2erase i];
    elseif ~strcmp(splitName{2} , typeHL) && ~strcmp(splitName{2} ,['2591' typeHL])
        i2erase = [i2erase i];
    end
    
end

% Erase non-wav file
fileList(i2erase)  = [];

% numer of file
nbF = length(fileList);

%Getting time from name
formatIn = 'yyyymmddThhMMss';
for i=1:nbF
    splitName = split(fileList{i}, '_');
    dateString = splitName{3};    
    dateN = datenum(dateString,formatIn);
    dateT(i,1) = datetime(dateN,'ConvertFrom', 'datenum');
    id{i} = splitName{4};
    
end

%Find the file
if ~exist('dateT')  
    error(['No file found for ' datestr(dateIn) ' in ' folderIn  '.'])
end


% Loading all time if not specified
if isempty(dateIn)
   dateIn =  dateT; 
end


for i=1:length(dateIn)
% Check if need to load another file 
bolTime =  dateT >= dateIn(i);
if any(bolTime) == 1
    ft = find(bolTime, 1, 'first');
    diffT = dateT(ft) - dateIn(i) ;
    if diffT == seconds(0)
        iFile(i) = ft;
    elseif diffT < minutes(5) && diffT > 0 && ft > 1
        iFile(i) = ft - 1 ;
    else
        error('Couldnt load the file at l.61 or add more conditions.')
    end
elseif dateT(end) + minutes(4.99) >  dateIn(i)
    iFile(i) = length(dateT);
else
    error('Couldnt load the file!')
end

% Check is file actually exist
%disp('Check function here in getWavName')
if iFile(i) == 0
    error('Error in bringRead.m. Request time is not in the folder.')
end
end

%iFile = unique(iFile);
%output name

if exist('iFile')
    fileName = fileList(iFile)';
    for i=1:length(iFile)
        splitName = split(fileList{iFile(i)}, '_');
        wi(i).arrLoc = splitName{1};
        wi(i).type = splitName{2};
        wi(i).dateStr = splitName{3};
        wi(i).wavID = splitName{4};
        
        dateString = splitName{3};    
        dateN = datenum(dateString,formatIn);
        wi(i).dateT = datetime(dateN,'ConvertFrom', 'datenum');
    end
else
        error('Error in bringRead.m. Request time is not in the folder.')
end
end




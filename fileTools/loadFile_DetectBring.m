% Load detection excel

% Path and file
path2File = ['/Users/Administrator/Documents/MPO/BRing/Data/detections/'];         
file2Load = ['BF_PRC_MLB_NARW_PR.xlsx'];

tab = readtable([path2File file2Load],'sheet','MLB');
nbF =  size(tab,1);

%tab.Properties.VariableNames{'FichierPNG'} = 'file';
%tab.Properties.VariableNames{'NARWTime_s_'} = 'time';

% Defautl parameter
formatIn = 'yyyymmddThhMMss';
spacer = '_';
%% Extract results

% Pre initialisation
for ii=1:nbF
   
    splitName = split(tab.file{ii}, '_');
    dateString = splitName{4};
    dateN = datenum(dateString,formatIn);
    dtime(ii) = datetime(dateN,'ConvertFrom', 'datenum') + seconds(tab.time(ii));
    
end

%%
save([path2File 'DETECT_NARW_MLB_PR_V1.mat'],'dtime')
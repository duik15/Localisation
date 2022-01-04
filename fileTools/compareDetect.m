% This scrip load detection and to make mapping
clear all
close all
path2Detect = '/Users/Administrator/Documents/MPO/BRing/Data/detections/';

ID{1} = 'PRC';
ID{2} = 'MLB';

ii=1;
file{ii} = ['DETECT_NARW_' ID{ii} '_PR_V1.mat'];
D(ii) = load([path2Detect file{ii}]);

ii=2
file{ii} = ['DETECT_NARW_' ID{ii} '_PR_V1.mat'];
D(2) = load([path2Detect file{2}]);

% Dafautl parameter
tDelayMax =  15;
% Modify structure
%save([path2Detect file{ii}],'-struct','D','dtime','angleAD','angleMD','energiedB','file')
%% Get the file name
for ii=1:2
    % Get file name
    D(ii).file = getWavName(D(ii).dtime, getDirectory(['toaster' ID{ii} ]))
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
%%


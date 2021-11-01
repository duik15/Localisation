function [index closeVal] = getClosest(matrix, wantedVal)
%Same as the min fonction of MATLAB but can take a vector and return only
%the index

%Preaalocaiton
index = nan(numel(wantedVal),1);

for i=1:numel(wantedVal)
[~,index(i)] = min(abs(matrix - wantedVal(i))); % Index of the depth 
end

index = unique(index);
closeVal = matrix(index);

end
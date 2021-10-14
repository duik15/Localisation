function [index matOut] = getInbetween(matIn, minMax)
% Get the index that fits between two values 
%[index matOut] = getInbetween(matIn, minMax)% 
%matIn = ptimeT
%minMax = [min(ptimeP) max(ptimeP)]

    ifirst = find(  0<(matIn - min(minMax)), 1 , 'first' ); 
    iend = find(  0<(matIn - max(minMax)), 1 , 'first' ); 
    
    if isempty(iend)
       iend = length(matIn); 
    end
    
    matOut = matIn(ifirst:iend);
    index = ifirst : iend;
    
end
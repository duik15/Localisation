function [loc offSet arrOri]= getArrLoc(hID)
% Get the most accurate array location for the moment
%[loc offSet arrOri]= getArrLoc(hID)
warning('getArrLoc.m is obsolet. Please use getArrInfo and recode your script using a structure.')

switch lower(hID)
    case 'aav'
        loc = [49.0907 -64.5372];
        offSet = -23;
        arrOri = 'clock';
    case 'cld'
        loc = [49.1933  -64.8315];
        offSet = 97;
        arrOri = 'counter';
    case 'mlb'
        loc = [48.6039 -64.1863];
        offSet = -78;
        arrOri = 'clock';
    case 'prc'
        loc= [48.5311 -64.1983];
        offSet = 135;
        arrOri = 'clock';
        
    otherwise
        error(['Can''t location. AV = Anse-a-Valleau, CD = Cloridorme, MB = Malbay, PC = Perce.' ])
end
        
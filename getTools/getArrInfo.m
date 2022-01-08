function [loc ai]= getArrInfo(hID)
% Get the most accurate array location for the moment
% [loc arr]= getArrInfo(hID)
% Where loc is the geopositon of the arr and arr is a structure containing
% many informations like the order, hydro. position,  
%[loc offSet arrOri]= getArrLoc(hID)  -> old Syntaxe obsolet.

switch lower(hID)
    case 'aav'
        ai.loc = [49.0907 -64.5372];
        ai.offSet = -23;
        ai.arrOri = 'clock';
        ai.xh = [0 2.8153 5.8708 8.8716 11.5519 13.4905 13.4905 14.5916 14.3387 13.6443  11.8551  9.4358  6.3276  3.7707  3.1079  0.5411 -2.5568 -4.9836 -4.9694 -2.4206];
        ai.yh = [0 0.8658 0.8435 -0.0394 -1.7619 -4.0724 -6.6749 -9.3099 -11.9984 -14.3325 -16.1202 -16.5980 -14.8928 -12.2318 -9.3739 -7.8348 -7.4210 -5.7031 -2.9183 -1.3566];
        ai.azimutMax = [-90, 90];
    case 'cld'
        ai.loc = [49.1933  -64.8315];
        ai.offSet = 97;
        ai.arrOri = 'counter';
        ai.xh=[0 -0.7903 0.7100 3.2029 6.0198 9.0990 10.9002 8.9040 6.2348 5.0490 6.4380 9.7425 12.2334 12.2306 10.1136 8.6284 8.7481 8.4264 5.8869 2.9986];
        ai.yh=[0 -2.7482 -5.384 -7.0614 -7.6449 -7.8366 -8.5376 -9.4106 -11.0053 -13.3869 -15.7872 -16.1861 -14.7605 -11.6937 -9.3713 -6.7706 -3.7894 -0.8292 1.1234 1.2759];
        ai.azimutMax = [-60, 120];
    case 'mlb'
        ai.loc = [48.6039 -64.1863];
        ai.offSet = -78;
        ai.arrOri = 'clock';
        ai.xh=[0 1.1940 3.0981 5.7060 8.6653 11.6901 14.5182 16.8067 18.4230 19.0958 18.9235 17.6643 15.7067 13.1221 10.1571 7.1448 4.3031 2.0291 0.4270 -0.2247];
        ai.yh=[0 2.7759 5.1765 6.8422 7.6537 7.5305 6.5924 4.7341 2.2917 -0.4581 -3.4326 -6.2230 -8.5628 -10.2413 -11.0769 -11.0141 -10.0422 -8.2079 -5.7896 -2.8977];
        ai.azimutMax = [0, 180];
    case 'prc'
        ai.loc= [48.5311 -64.1983];
        ai.offSet = 135;
        ai.arrOri = 'clock';
        ai.xh=[0 -2.0409 -4.6008 -7.4667 -10.2984 -12.8242 -14.7966 -15.9877 -16.2946 -15.5958 -14.7413 -12.1864 -9.6619 -8.1987 -6.7088 -3.9947 -1.1311 1.6361 3.0529 1.5787]; 
        ai.yh=[0 -1.1660 -2.5307 -2.5352 -1.5951 0.1105 2.5024 5.3343 8.2865 10.9606 13.5860 14.5295 12.9216 10.1355 7.5359 6.2587 6.4663 7.1628 4.5436 2.1991];
        ai.azimutMax = [0, 180];
    otherwise
        error(['Can''t location. AV = Anse-a-Valleau, CD = Cloridorme, MB = Malbay, PC = Perce.' ])
end % end switch

% First outarg: Will change eventually with updated code.
loc = ai.loc;

% Center arr to center after new position have been input.
x0 = mean(ai.xh);
y0 = mean(ai.yh);
ai.xh = ai.xh - x0;
ai.yh = ai.yh - y0;

% Sound Speed 
ai.c = 1475;

end 
        
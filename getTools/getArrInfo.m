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
        ai.xh = [0 2.7033 5.6089 8.4344 10.5549 12.3247 13.2599 13.3094 12.3194 10.8610  8.8097  6.1495  4.1051  3.4854  2.9988  0.6787 -2.2052 -4.7213 -4.8024 -2.2494];
        ai.yh = [0 0.98115 1.0462 0.23029 -1.1796 -3.1658 -5.7182 -8.3689 -11.286 -13.5013 -15.4163 -15.9193 -14.4628 -11.4901 -9.0871 -7.6885 -7.509 -6.0617 -3.7395 -1.6627];
        ai.azimutMax = [-90, 90];
    case 'cld'
        ai.loc = [49.1933  -64.8315];
        ai.offSet = 97;
        ai.arrOri = 'counter';
        ai.xh=[0 -0.85442 0.50753 2.9847 5.7683 8.7395 10.3849 8.2085 6.0404 4.831 6.0476 8.9899 11.413 11.199 9.3466 7.8237 7.7304 7.2501 5.2569 2.5475];
        ai.yh=[0 -2.7482 -5.384 -7.0614 -7.6449 -7.8366 -8.5376 -9.4106 -11.0053 -13.3869 -15.7872 -16.1861 -14.7605 -11.6937 -9.3713 -6.7706 -3.7894 -0.8292 1.1234 1.2759];
        ai.azimutMax = [-60, 120];
    case 'mlb'
        ai.loc = [48.6039 -64.1863];
        ai.offSet = -78;
        ai.arrOri = 'clock';
        ai.xh=[0 0.91706 2.4702 4.707 7.3827 10.3162 13.159 15.5768 17.5015 18.2361 18.4047 17.7424 16.1233 13.7894 11.0051 8.0391 5.1927 2.7985 0.98695 0.076049];
        ai.yh=[0 2.8335 5.2546 6.8999 7.8402 7.8818 7.1528 5.473 3.1381 0.19332 -2.6932 -5.5842 -7.9302 -9.5627 -10.4316 -10.4667 -9.6166 -7.9123 -5.718 -2.9604];
        ai.azimutMax = [0, 180];
    case 'prc'
        ai.loc= [48.5311 -64.1983];
        ai.offSet = 135;
        ai.arrOri = 'clock';
        ai.xh=[0 -2.2134 -5.0383 -8.0563 -10.9186 -13.2934 -15.3135 -16.2057 -15.9951 -15.1329 -14.1683 -11.6308 -9.0821 -7.5243 -6.2387 -3.6652 -0.83313 1.9931 3.3276 1.5583];
        ai.yh=[0 -1.9489 -3.2505 -3.5494 -2.8135 -1.2863 0.86294 3.6998 6.6637 9.484 12.157 13.5441 12.1563 9.2514 6.8973 5.9277 6.4666 7.3163 4.9562 2.2587];
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
        
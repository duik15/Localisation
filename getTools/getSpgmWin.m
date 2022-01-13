function spgm = getSpgmWin(spgm,varargin)
   
%%
if ~isempty(varargin)
        audioInfo = varargin{1};   
        % Use audioInfo to get more information
        spgm.fs = audioInfo.Fs; spgm.ns =  audioInfo.Ns; spgm.im.ns= spgm.ns; spgm.im.dur = audioInfo.dura; 
end

% Spectrogram window parameter
    spgm.win.ns = fix(spgm.win.dur*spgm.fs);
    if strcmp(spgm.win.type , 'hanning')
        spgm.win.val = hanning(win.ns);
    elseif strcmp(spgm.win.type , 'kaiser')
        spgm.win.wpond = kaiser(spgm.win.ns ,0.1102*(180-8.7));
        spgm.win.val = spgm.win.wpond*sqrt(spgm.win.ns/sum(spgm.win.wpond.^2)); %w_pond
    end
    spgm.win.novlp = fix(spgm.win.ovlp/100*spgm.win.ns);
    spgm.win.nfft = max(256,2^nextpow2(spgm.win.ns*10));
    
    spgm.im.fmin = spgm.im.freqlims(1);
    spgm.im.fmax = spgm.im.freqlims(2);
    
    
end
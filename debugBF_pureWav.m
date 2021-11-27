%% Variables
% enviroment variables
arrID = 'AAV';
mainCh = [4 10 15];    % Looking only one chanel to speed up everything
%outName = 'aavLiveBoat_onlyTrackTime_Ns14_f150-200hz';
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderOut = [folderIn '/' 'spectrogram/'];
%folderIn = ['E:\Bring_Dep_1\' arrID '\'];
folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\debug\']; %spectrogram'];

saveflag = true;                   % [true or false]


% Selecting file wanted
time2Load = datetime(2021,07,15,18,05,00);
%time2Load = datetime(2021,07,14,00,00,00):minutes(5):datetime(2021,07,14,02,00,00);

% List file in the folderIn
file = getWavName(time2Load,folderIn);


% spectrogram image parameters
spgm.im.freqlims = [10 200];       % [Hz] frequency scale boundary limits
spgm.im.freqscale = 'linear';       % ['linear' or 'log'] frequency scale type
spgm.im.clims = [-120 0];           % [dB] C limite pcolor
spgm.im.lbwh = 'auto';              % [left, bottom, width, height or 'auto'] image size and location [800 500 188 179];
spgm.im.format = 'png';             % image output format 'png', 'jpg'
spgm.im.dur = 30;%'all';                % [s or 'all'] figure duration
spgm.im.ovlp = 0;                   % [%] image window overlap
spgm.im.link_axes = false;           % [true false] use fonction linkaxes
spgm.im.figvision = true  ;           % [true false] visiblity of figure before saveas jpf file
spgm.im.constrast = [1 1.2];        % Constrat the clims by [x x] between 0 and 1 caxis = caxis * [cont1 cont2]
spgm.im.movm = 100;                  % Movmean parameter for the first panel
spgm.im.FontS = 14;                  % Image font Size
% spectrogram window parameters
spgm.win.dur = 2048/10000;                  % [s] spectrogram window length duration
spgm.win.ovlp = 50;                 % [%] spectrogram window overlap
spgm.win.type = 'hanning';
spgm.win.opad = 8;                  % [s] spectrogram zero padding


if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

% Figure parameter
sp.height=15; sp.width=35; sp.nbx=1; sp.nby= length(mainCh)+1;
sp.ledge=4; sp.redge=6; sp.tedge=1.5; sp.bedge=3;
sp.spacex=0.5; sp.spacey=0.5;
sp.fracy = [0.1 0.3 0.3 0.3];
sp.pos=subplot2(sp);




%% Open the first
acinfo = audioinfo([folderIn file{1}]);
[acf.s,acf.fs] = audioread([folderIn file{1}]);


acf.s = acf.s(:,mainCh);
if strcmp(spgm.im.dur,'all')
    spgm.im.dur = acinfo.Duration;
end

spgm.im.ns = fix(spgm.im.dur*acinfo.SampleRate);
spgm.im.n = round(acinfo.TotalSamples/spgm.im.ns);

% define spectrogram window
win.ns = fix(spgm.win.dur*acf.fs);
win.val = hanning(win.ns);
win.novlp = fix(spgm.win.ovlp/100*win.ns);
win.nfft = fix((spgm.win.dur+spgm.win.opad)*acf.fs);


%%
%acinfo = audioinfo([folderIn file{1}]);
%[acf.s,acf.fs] = audioread([folderIn file{1}]);
acinfo = audioinfo([folderIn file{1}]);
winL = spec.dur * acinfo.SampleRate;
[wav.s,wav.fs] = audioread([folderIn file{1}],[1 winL]);
spec.Fs = wav.fs;
spec.Ns = size(wav.s,1);
acf = wav;

mainCh = 1;
acf.s = acf.s(:,mainCh);
if strcmp(spgm.im.dur,'all')
    spgm.im.dur = acinfo.Duration;
end

spgm.im.ns = fix(spgm.im.dur*acinfo.SampleRate);
spgm.im.n = round(acinfo.TotalSamples/spgm.im.ns);

% define spectrogram window
win.ns = fix(spgm.win.dur*acf.fs);
win.val = hanning(win.ns);
win.novlp = fix(spgm.win.ovlp/100*win.ns);
win.nfft = fix((spgm.win.dur+spgm.win.opad)*acf.fs);

i_file=1
acinfo = audioinfo([folderIn file{i_file}]);
    [acf.s,acf.fs] = audioread([folderIn file{i_file}]);
    
    j=1
    
    
    
    disp(['Inter no: ' num2str(j)])
        if j == spgm.im.n
            (j-1)*spgm.im.dur * acinfo.SampleRate+1 : j*spgm.im.dur * acinfo.SampleRate;
        else
            it = (j-1)*spgm.im.dur * acinfo.SampleRate+1 : j*spgm.im.dur * acinfo.SampleRate;
        end
        %
        timeRec = it / acinfo.SampleRate;
        
        % HERE
        
        %acf.s = acf.s(:,mainCh);
        rec = acf.s(it,mainCh);
        
        % make spectrogram
        disp([datestr(datetime('now')) ' | Make spectrogram ']);
        i_ch=1%:length(mainCh)
            val = win.val;
            novlp = win.novlp;
            mfft=win.nfft;
            %[~,F,T,P] = spectrogram(rec(:,i_ch),win.val,win.novlp,win.nfft,acf.fs);
            [~,F,T,P] = spectrogram(matWav(:,1),win.val,win.novlp,win.nfft,spec.Fs); ok1=1
            PdB(:,:,i_ch) = 10*log10(P);
        
        timeP = spgm.im.dur * (j-1) + T;
        %%
        
        FF = freqV;
        TT = timeV;
        MM = mat1;
        
        %%
        
        % Panel one
        %ax(1)=axes('position',sp.pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',spgm.im.FontS,'Box','on','Layer','top');
        
        %plot(timeRec,movmean(rec,spgm.im.movm,1))
        %ax(1).XTick = [];
        %ylabel('Amp.')
        %datetick('x','hh:mm:ss')
        clf
        figure(1)
        % Panel two
        % pcolor is slightly slower than the surf function but simpler to use
        %for i_ch=1:length(mainCh)
         i_ch=1   
            pcolor(TT,FF,MM)
            %pcolor(timeV,freqV,mat1)
            %pcolor(timeP,F,PdB(:,:,i_ch));
            shading flat;
            axis tight;
            ylim(spgm.im.freqlims);
            set(gca,'YScale',spgm.im.freqscale);
            %caxis([-120 round(max(PdB(:)),-1)+10] .* spgm.im.constrast)
            colormap(jet)
            %set(gca,'Layer','top')
            grid on
            set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
            % Text and label
            text(timeP(end) - 3 , spgm.im.freqlims(2) * 0.9,['ch.' num2str(mainCh(i_ch))], 'backgroundcolor', [1 1 1],'margin',2)
            %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2 + 5,['ch.'])
            %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2-5,[num2str(mainCh(i_ch))])
            %if i_ch == ceil(sp.nby/2)
                yl = ylabel('Frequency (Hz)','FontSize',spgm.im.FontS);
            %end
            
        %end
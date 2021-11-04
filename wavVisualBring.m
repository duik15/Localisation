%% WAV to fig spectrogram
%
% DESCRIPTION:
% Load audio file, make spectrogram and convert fig in jpg file
% Use IML-PAM-Toolbox\wav_tools\wav_spectrogram_image.m to make this script developed by Florian A and Clement J
%
%
% UPDATES:
% 2021-09-08:       Clement.Juif@dfo-mpo.gc.ca (CJ)
% 2021-11-03        Kevin.duquette@dfo-mpo.gc.ca (KD)
%
%

clc
clear all
close all
tic
%% Variables
% enviroment variables
arrID = 'AAV'; 
mainCh = 2;    % Looking only one chanel to speed up everything
%outName = 'aavLiveBoat_onlyTrackTime_Ns14_f150-200hz';
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
folderOut = [folderIn '/' 'spectrogram/'];
saveflag = true;                   % [true or false]


% Selecting file wanted
time2Load = datetime(2021,07,15,18,10,00):minutes(5):datetime(2021,07,15,18,25,00);

% List file in the folderIn
file = getWavName(time2Load,folderIn);


% spectrogram image parameters
spgm.im.freqlims = [10 500];       % [Hz] frequency scale boundary limits
spgm.im.freqscale = 'linear';       % ['linear' or 'log'] frequency scale type
spgm.im.clims = [-120 0];           % [dB] C limite pcolor
spgm.im.lbwh = 'auto';              % [left, bottom, width, height or 'auto'] image size and location [800 500 188 179];
spgm.im.format = 'png';             % image output format 'png', 'jpg'
spgm.im.dur = 'all';                % [s or 'all'] figure duration
spgm.im.ovlp = 0;                   % [%] image window overlap
spgm.im.link_axes = false;           % [true false] use fonction linkaxes
spgm.im.figvision = true;           % [true false] visiblity of figure before saveas jpf file

% spectrogram window parameters
spgm.win.dur = 0.1;                  % [s] spectrogram window length duration
spgm.win.ovlp = 50;                 % [%] spectrogram window overlap
spgm.win.type = 'hanning';
spgm.win.opad = 0;                  % [s] spectrogram zero padding


if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

% Figure parameter
sp.height=20; sp.width=30; sp.nbx=1; sp.nby=2;
sp.ledge=3; sp.redge=5; sp.tedge=1.5; sp.bedge=3;
sp.spacex=0.5; sp.spacey=0.5; sp.fontsize=16;
sp.fracy = [0.2 0.8];
sp.pos=subplot2(sp);




%% Open the first
acinfo = audioinfo([folderIn file{1}]);
[acf.s,acf.fs] = audioread([folderIn file{1}]);

timeFile = getFileTime([folderIn file{1}]);

% Create a time vector
time =  timeFile + seconds((1:acinfo.TotalSamples)/acinfo.SampleRate ) ;

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

%% process
for i_file = 1:length(file)
    % Read file
    disp([datestr(datetime('now')) ' | Load file -> ' file{i_file}]);
    acinfo = audioinfo([folderIn file{i_file}]);
    [acf.s,acf.fs] = audioread([folderIn file{i_file}]);
    
    acf.s = acf.s(:,mainCh);
    
    % make spectrogram
    disp([datestr(datetime('now')) ' | Make spectrogram ']);
    if spgm.im.figvision
        fig_hdl = figure();
    else
        fig_hdl = figure('visible','off');
    end
    
    i_ch=1; % Skip a foorlop on chanel
        % classical spectrogram
        % t1_spectro = tic;
        [~,F,T,P] = spectrogram(acf.s(:,i_ch),win.val,win.novlp,win.nfft,acf.fs);
        PdB = 10*log10(P);
        Tdt = datetime(2021,07,15,18,10,00) + seconds(T);
        % t2_spectro = toc(t1_spectro);
        
        % t1_plotfig = tic;
        %ax(i_ch) = subplot(size(acf.s,2)/2,2,i_ch);
        
        % --------------- Figure --------------------------
        set(fig_hdl, 'PaperUnits', 'centimeters');
        set(fig_hdl, 'PaperSize', [sp.width sp.height]);
        set(fig_hdl, 'PaperPositionMode', 'manual');
        %ax(i_ch) = subplot(size(acf.s,2),2,i_ch);
        
        % Panel one
        ax(1)=axes('position',sp.pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',sp.fontsize,'Box','on','Layer','top');
        plot(time,acf.s)
        ax(1).XTick = [];
        ylabel('Amplitude')
        %datetick('x','hh:mm:ss')
        
        % Panel two
        % pcolor is slightly slower than the surf function but simpler to use
        ax(2)=axes('position',sp.pos{1,2},'XGrid','off','XMinorGrid','off','FontSize',sp.fontsize,'Box','on','Layer','top');
        pcolor(Tdt,F,PdB);
        shading flat;
        axis tight;
        ylim(spgm.im.freqlims);
        set(gca,'YScale',spgm.im.freqscale);
        %caxis(spgm.im.clims);
        caxis([-120 round(max(PdB(:)),-1)+10])
        colormap(jet)
        ylabel(num2str(i_ch))
        grid on
        set(gca,'Layer','top')
        ylabel('Frequecy (Hz)')
        
        % Color bar
        cbPos = [sp.pos{2}(1)+sp.pos{2}(3)+sp.spacex sp.pos{2}(2) 0.03 sp.pos{2}(4)+sp.pos{1}(4)+sp.spacey ];
        cb= colorbar;%('Position',cbPos);
        %cb.Position = cbPos
        ax(2).Position = sp.pos{2};
        ylabel(cb,'Power (dB)')
        
        if ~strcmp(spgm.im.lbwh,'auto')
            set(fig_hdl,'Position',spgm.im.lbwh);
        end
        
        sgtitle(file{i_file},'Interpreter','none');
        
        if spgm.im.link_axes
            linkaxes(ax,'xy');
        end
        
        if saveflag
            disp([datestr(datetime('now')) ' | save file ' file{i_file}(1:end-3) spgm.im.format ' in the folder ' folderOut]);
            saveas(fig_hdl,[folderOut 'test11' file{i_file}(1:end-3) spgm.im.format]);
        end
        if spgm.im.figvision
        pause;
        end
        
        clf
        %clear acf fig_hdl
    end % .. i_file
toc

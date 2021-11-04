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
arrID = 'AAV'; mainCh = 2;
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
spgm.im.format = 'jpg';             % image output format 'png', 'jpg'
spgm.im.dur = 'all';                % [s or 'all'] figure duration
spgm.im.ovlp = 0;                   % [%] image window overlap
spgm.im.link_axes = false;           % [true false] use fonction linkaxes
spgm.im.figvision = false;           % [true false] visiblity of figure before saveas jpf file

% spectrogram window parameters
spgm.win.dur = 5;                  % [s] spectrogram window length duration
spgm.win.ovlp = 50;                 % [%] spectrogram window overlap
spgm.win.type = 'hanning';
spgm.win.opad = 0;                  % [s] spectrogram zero padding


if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

% Figure parameter
sp.height=20; sp.width=30; sp.nbx=1; sp.nby=2;   
sp.ledge=2; sp.redge=3; sp.tedge=1.2; sp.bedge=2;
sp.spacex=0.5; sp.spacey=0.5; sp.fontsize=16; 

%% process
for i_file = 1%:length(file)
    %% Read file
    disp([datestr(datetime('now')) ' | Load file -> ' file{i_file}]);
    acinfo = audioinfo([folderIn file{i_file}]);
    [acf.s,acf.fs] = audioread([folderIn file{i_file}]);
    
    acf.s = acf.s(:,mainCh);
    if strcmp(spgm.im.dur,'all')
        spgm.im.dur = acinfo.Duration;
    end
    
    spgm.im.ns = fix(spgm.im.dur*acinfo.SampleRate);
    spgm.im.n = round(acinfo.TotalSamples/spgm.im.ns);
    
    %% define spectrogram window
    win.ns = fix(spgm.win.dur*acf.fs);
    win.val = hanning(win.ns);
    win.novlp = fix(spgm.win.ovlp/100*win.ns);
    win.nfft = fix((spgm.win.dur+spgm.win.opad)*acf.fs);
    
    %% make spectrogram
    disp([datestr(datetime('now')) ' | Make spectrogram ']);
    if spgm.im.figvision
        fig_hdl = figure();
    else
        fig_hdl = figure('visible','off');
    end
    
    for i_ch = 1:size(acf.s,2) 
        % classical spectrogram
        % t1_spectro = tic;
        [~,F,T,P] = spectrogram(acf.s(:,i_ch),win.val,win.novlp,win.nfft,acf.fs);
        PdB = 10*log10(P);
        % t2_spectro = toc(t1_spectro);
        
         % t1_plotfig = tic;
        %ax(i_ch) = subplot(size(acf.s,2)/2,2,i_ch);
        
        % --------------- Figure --------------------------
        %f=figure(1);
        %clf(f);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [sp.width sp.height]);
        set(gcf, 'PaperPositionMode', 'manual');
        %ax(i_ch) = subplot(size(acf.s,2),2,i_ch);
        
        % DEV NOTE: 
        % pcolor is slightly slower than the surf function but simpler to use
        ax_hdl = pcolor(T,F,PdB);
        shading flat;
        axis tight;
        ylim(spgm.im.freqlims);
        set(gca,'YScale',spgm.im.freqscale);
        set(ax_hdl,'lineStyle','none');
        if i_ch ~= size(acf.s,2) && i_ch ~= size(acf.s,2)-1
            set(gca,'XTick',[]);
        else
            xlabel('Time [s]')
        end
        caxis(spgm.im.clims);
        colormap(jet)
        ylabel(num2str(i_ch))
        grid on
        set(gca,'Layer','top')
%         t2_plotfig = toc(t1_plotfig);
        
    end % .. i_ch
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
    
    close all
    clear acf fig_hdl
end % .. i_file
t2=toc

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
infoGL = opengl('info');
opengl('software')
% enviroment variables
arrID = 'PRC';
mainCh = [7 3 19];    % Looking only one chanel to speed up everything
%outName = 'aavLiveBoat_onlyTrackTime_Ns14_f150-200hz';
%folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderOut = [folderIn '/' 'spectrogram/'];
%folderIn = ['E:\Bring_Dep_1\' arrID '\'];
%folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\spectro\']; %spectrogram'];
folderIn = ['\\169.254.47.215\usbshare1-2\Bring_Dep_2\' arrID '\']
folderOut = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\' arrID '\spectro\']; %spectrogram'];

saveflag = true;                   % [true or false]


% Selecting file wanted
%time2Load = datetime(2021,07,15,09,50,00);
time2Load = datetime(2021,08,04,00,00,00):minutes(5):datetime(2021,08,05,00,00,00);

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
spgm.im.figvision = false  ;           % [true false] visiblity of figure before saveas jpf file
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

%% process
tic
for i_file = 1:length(file)
    % Read file
    %disp([datestr(datetime('now')) ' | Load file ' num2str(i_file) '/' num2str(length(file)) ' -> ' file{i_file}]);
    disp( [datestr(datetime('now')) ' | Treating spectrograme for file ' num2str(i_file) ' / ' num2str(length(file))])
    acinfo = audioinfo([folderIn file{i_file}]);
    [acf.s,acf.fs] = audioread([folderIn file{i_file}]);
    
    % Create a time vector
    %timeFile = getFileTime([file{i_file}]);
    %time =  timeFile + seconds((1:acinfo.TotalSamples)/acinfo.SampleRate ) ;
    %time = timeFile + seconds(it);
    
    % Loop on short
    for j=1:spgm.im.n %j=1:spgm.im.n
        %disp(['Inter no: ' num2str(j)])
        if j == spgm.im.n
            (j-1)*spgm.im.dur * acinfo.SampleRate+1 : j*spgm.im.dur * acinfo.SampleRate;
        else
            it = (j-1)*spgm.im.dur * acinfo.SampleRate+1 : j*spgm.im.dur * acinfo.SampleRate;
        end
        %
        timeRec = it / acinfo.SampleRate;
        
        %acf.s = acf.s(:,mainCh);
        rec = acf.s(it,mainCh);
        
        % make spectrogram
        %disp([datestr(datetime('now')) ' | Make spectrogram ']);
        for i_ch=1:length(mainCh)
            [~,F,T,P] = spectrogram(rec(:,i_ch),win.val,win.novlp,win.nfft,acf.fs);
            PdB(:,:,i_ch) = 10*log10(P);
        end
        timeP = spgm.im.dur * (j-1) + T;
        % A AJOUTER:
        % Filtre de bande passante
        %ordfilt=3; %ordre filt. butterworth
        %fbnd=[fmin,fmax];    %rmq: si iflim, fmin=flim(1), etc.; ou presque. Donc, ceci valide pour ifmin=0 ou 1
        %[Bbut,Abut]=butter(ordfilt,fbnd/(fsaud/2)); % Wn:  the normalized cutoff frequency Wn must be a number between 0 and 1,
        %      where 1 corresponds to the Nyquist frequency(fs/2), π radians per sample.
        %Wavaud = filter(Bbut,Abut,Wavaud);
        
        
        % Make fig
        if spgm.im.figvision
            fig_hdl = figure(1);
        else
            fig_hdl = figure('visible','off');
        end
        cmapBW = flipud(gray(128));
        
                
        % --------------- Figure --------------------------
        set(fig_hdl, 'PaperUnits', 'centimeters');
        set(fig_hdl, 'PaperSize', [sp.width sp.height]);
        set(fig_hdl, 'PaperPositionMode', 'manual');
        %ax(i_ch) = subplot(size(acf.s,2),2,i_ch);
        
        % Panel one
        ax(1)=axes('position',sp.pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',spgm.im.FontS,'Box','on','Layer','top');
        plot(timeRec,movmean(rec,spgm.im.movm,1))
        ax(1).XTick = [];
        ylabel('Amp.')
        %datetick('x','hh:mm:ss')
        
        % Panel two
        % pcolor is slightly slower than the surf function but simpler to use
        for i_ch=1:length(mainCh)
            ax(i_ch+1)=axes('position',sp.pos{1,i_ch+1},'XGrid','off','XMinorGrid','off','FontSize',spgm.im.FontS,'Box','on','Layer','top');
            pcolor(timeP,F,PdB(:,:,i_ch));
            shading flat;
            axis tight;
            ylim(spgm.im.freqlims);
            set(gca,'YScale',spgm.im.freqscale);
            caxis([-120 round(max(PdB(:)),-1)+10] .* spgm.im.constrast)
            colormap(cmapBW)
            %set(gca,'Layer','top')
            grid on
            set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
            % Text and label
            text(timeP(end) - 3 , spgm.im.freqlims(2) * 0.9,['ch.' num2str(mainCh(i_ch))], 'backgroundcolor', [1 1 1],'margin',2)
            %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2 + 5,['ch.'])
            %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2-5,[num2str(mainCh(i_ch))])
            if i_ch == ceil(sp.nby/2)
                yl = ylabel('Frequency (Hz)','FontSize',spgm.im.FontS);
            end
            
        end
        
        
        % Color bar
        cbPos = [sp.pos{end}(1)+sp.pos{end}(3)+0.01 sp.pos{end}(2) 0.03 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
        cb= colorbar;%('Position',cbPos);
        cb.Position = cbPos;
        ax(end).Position = sp.pos{end};
        ylabel(cb,'Power (dB)')
        
        % General text and label
        sgtitle(file{i_file},'Interpreter','none');
        xlabel('Time related to file (s)','FontSize',spgm.im.FontS)
        
        if ~strcmp(spgm.im.lbwh,'auto')
            set(fig_hdl,'Position',spgm.im.lbwh);
        end
        
        
        if spgm.im.link_axes
            linkaxes(ax,'xy');
        end
     
        if saveflag
            %disp([datestr(datetime('now')) ' | save file ' file{i_file}(1:end-4) '_s' num2str(floor(timeRec(1))) '_e' num2str(ceil(timeRec(end))) ' in the folder ' folderOut]);
            print('-dpng','-r150',[folderOut file{i_file}(1:end-4) '_s' num2str(floor(timeRec(1))) '_e' num2str(ceil(timeRec(end))) '.' spgm.im.format])
            %saveas(fig_hdl,[folderOut file{i_file}(1:end-4) '_s' num2str(floor(timeRec(1))) '_e' num2str(ceil(timeRec(end))) '.' spgm.im.format]);
        end
           
        if spgm.im.figvision
            pause;
            clf
        end


    end
      disp('End file')
      toc
    %clear acf fig_hdl
end % .. i_file
disp('End of script')
toc

infoGL

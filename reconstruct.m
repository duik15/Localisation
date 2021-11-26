function [reconFFT, timeV, freqVFr ] = reconstruct(arrIDorMATPOND, matWav,azimut, spec, varargin)
%[recon, timeV, freqVFr ] = reconstruct(nbV,matPondahf, matFFT, azimut, freq, spec)
% This fonction reconstruct signal il all different direction using all
% chanel. So all signal reconstruct point in one direction of listening.
% Input var: matPondahf = Ponderation matrix with indice (nbAzimut,nbMicro, nbFreq)
%            matFFT     = matric of spectrogram (nbFreq, nbTimefft)
%            azimut = vecteur of azimut
%            spec =  sepctogram parameter
% Debug:
%matPondahf = MAT_POND_vs_azim_h_freq;
%matWav = wav.sdb;
%azimutV = (0:nbV-1) * 360/nbV;
%nbV = length(azimutV);

%%
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'fmin'
                fmin = varargin{2};
            case 'fmax'
                fmax = varargin{2};
            case 'matpond'
                matPondahf = varargin{2};
            case 'arrid'
                arrID = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end
%% Defaut value

if ischar(arrIDorMATPOND)
   arrID =  arrIDorMATPOND;
else
    matPondahf = arrIDorMATPOND;
end

% Default parameter
showFig = [0];
Ns = size(matWav,1);
nbV = length(azimut);


% Matrice of ponderation
if ~exist('matPondahf')
    % Get the ponderation matric
    matPondahf = makePond(arrID,azimut,Ns,spec.Fs,'fmin',spec.fmin,'fmax',spec.fmax);
    %matPondahf = makePond(arrID,1:360,Ns,spec.Fs,'fmin',spec.fmin,'fmax',spec.fmax);
end
%% Create the fft of the signal

% Create the frequence vecteur relate to FFT
freq= (0:1:Ns-1)*spec.Fs/Ns;
indF = find( (freq >= spec.fmin)&(freq <= spec.fmax));
freqF = freq(indF);

% Create matFFT with only freq needed
matFFT2 = nan(size(matWav));
for ii = 1:size(matWav,2)
    matFFT2(:,ii) = fft(matWav(:,ii));
end
matFFT = matFFT2(indF,:); clear matFFT2;


%% Recontruct wav signal and its FFT

% Initialize
recon=zeros(spec.Ns,nbV);
for u = 1 : nbV
    % Initiate
    matFFTPond = zeros(1,floor(Ns/2)+1);
    % Calcul the fft for the direction and add those freq values to the null matrix
    matFFTPondSum = ( sum(matFFT .* squeeze(matPondahf(u,:,:))',2 ))';
    
    % Add the conj
    matFFTPond(indF) = matFFTPondSum;
    matFFTPond2 = [ matFFTPond conj(fliplr(matFFTPond(2:length(matFFTPond)-1))) ];
    %matFFTPondV(:,u) =  matFFTPond2';
    % Reconstruct WAV % This matrice could be remove to save memory and
    % time
    recon(:,u) = ifft(matFFTPond2);
    
    % Calcul of spectrogrammes associated with BF
    [timeV, freqV, ~ , matFFTV_dB] = COMP_STFT_snapshot(recon(:,u),0, spec.Fs, spec.winSz, spec.rec, spec.wpond, spec.zp);
    indFr = find( ((freqV>= spec.fmin)&(freqV <= spec.fmax)));
    reconFFT(u,:,:)=matFFTV_dB(:,indFr);
end


freqVFr = freqV(indFr);

%%
if showFig ==1
    % ------ Figure plot --------
    cmapBW = flipud(gray(128));
    figure(1)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [sp.width sp.height]);
    set(gcf, 'PaperPositionMode', 'manual');
    
    meshV = reshape(1:nbV,sp.nbx,sp.nby);
    
    for u = 1 :  nbV
        [sx,sy] = find(meshV==u);
        ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
        %subplot(nbSub,nbSub,u)
        M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
        pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat;
        ylim([fmin_int fmax_int])
        
        if sx==1
            ylabel(' f (Hz)')
        else
            ax(u).YTickLabel = [];
        end
        if sy==sp.nby
            xlabel(' t (s)')
        else
            ax(u).XTickLabel = []
        end
        
        title(['V ' num2str(u) 'Az (°) ' num2str(floor(180/pi*vec_azim_nvoie(u)))])
        caxis([Lmin Lmax])
        colormap jet
        %colormap(cmapBW)
        grid on
        set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
    end
    
    sgtitle(fileList{iFile})
    
    
    if printFig ==true
        print([folderOut 'subplotSpectroNVoie_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f5')
    end
    
end


end

function [reconFFT, timeV, freqV, recon, matFFT] = beamForming(arrIDorMATPOND, matWav,azimut, spgm ,varargin)
%[recon, timeV, freqVFr ] = reconstruct(nbV,matPondahf, matFFT, azimut, freq, spec)
% This fonction reconstruct signal il all different direction using all
% chanel. So all signal reconstruct point in one direction of listening.
% Input var: matPondahf = Ponderation matrix with indice (nbAzimut,nbMicro, nbFreq)
%            matFFT     = matric of spectrogram (nbFreq, nbTimefft)
%            azimut = vecteur of azimut
%            spgm =  sepctogram parameter
% WARNING : THERE A BUG IN THE WAV RECONSTRUCTION. NEED TO BE FIX.
% Debug:
%matWav = wav.sdb;
%azimut = (0:nbV-1) * 360/nbV;
%nbV = length(azimutV);
specMethod  = 'spectro';
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
            case 'specmethod'
                specMethod = varargin{2};            
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
if ~exist('matPondahf') || size(matPondahf,2) ~= spgm.im.ns
    % Get the ponderation matric
    matPondahf = makePond(arrID,azimut,spgm.im.ns,spgm.fs,'fmin',spgm.im.fmin,'fmax',spgm.im.fmax);
    %matPondahf = makePond(arrID,1:360,Ns,spec.Fs,'fmin',spec.fmin,'fmax',spec.fmax);
end
%% Create the fft of the signal

% Create the frequence vecteur relate to FFT
freq= (0:1:spgm.im.ns-1)*spgm.fs/spgm.im.ns;
indF = find( (freq >= spgm.im.fmin)&(freq <= spgm.im.fmax));
freqF = freq(indF);

% Create matFFT with only freq needed
matFFT2 = nan(size(matWav));
for ii = 1:size(matWav,2)
    matFFT2(:,ii) = fft(matWav(:,ii));
end
matFFT = matFFT2(indF,:); clear matFFT2;


%% Recontruct wav signal and its FFT

% Initialize

if mod(spgm.im.ns,2) ==1
spgm.im.ns = spgm.im.ns-1;
matWav = matWav(1:end-1,:); 
end

recon=zeros(spgm.im.ns,nbV);
for u = 1 : nbV
    % Initiate
    matFFTPond = zeros(1,floor(Ns/2)+1);
    % Calcul the fft for the direction and add those freq values to the null matrix
    matFFTPondSum = ( sum(matFFT .* squeeze(matPondahf(u,:,:))',2 ))';
    matFFTPond(indF) = matFFTPondSum;
    % Add the conj
    matFFTPond2 = [ matFFTPond conj(fliplr(matFFTPond(2:length(matFFTPond)-1))) ];
    
    % Reconstruct WAV % This matrice could be remove to save memory and
    %size(matFFTPond2)
    %size(recon)
    %size(ifft(matFFTPond2))
    recon(:,u) = ifft(matFFTPond2);
    
    % Calcul of spectrogrammes associated with BF
    if strcmp(specMethod, 'spectro')
        [~,freqFV,timeV,tmp] = spectrogram(recon(:,u),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
        matFFTV_dB = 10*log10(tmp);
        
        % Keep only wanted frequency
        indFV = find( ((freqFV>= spgm.im.fmin)&(freqFV <= spgm.im.fmax)));
        reconFFT(u,:,:)=matFFTV_dB(indFV,:)';
    elseif strcmp(specMethod ,'cedric')
        disp('Using cedric method for spectrogram')
        % Calcul of spectrogrammes associated with BF
        [timeV, freqFV, ~ , matFFTV_dB] = COMP_STFT_snapshot(recon(:,u),0, spgm.fs, spgm.win.ns, spgm.win.ovlp/100, spgm.win.val, spgm.win.opad);
        indFV = find( ((freqFV>= spgm.im.fmin)&(freqFV <= spgm.im.fmax)));
        reconFFT(u,:,:)=matFFTV_dB(:,indFV);
    else
        error('Do not understant specMethod property.')
    end
end

freqV = freqFV(indFV);

% Flip the matrix because it need to be
reconFFT = flipdim(reconFFT,2);


%%
%{
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
        
        title(['V ' num2str(u) 'Az (?) ' num2str(floor(180/pi*vec_azim_nvoie(u)))])
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

%}
end

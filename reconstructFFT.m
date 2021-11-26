function [reconFFT, timeV, freqVFr ] = reconstructFFT(nbV,matPondahf, matFFT, azimut, freq, spec)
%[recon, timeV, freqVFr ] = reconstruct(nbV,matPondahf, matFFT, azimut, freq, spec)
% This fonction reconstruct signal il all different direction using all
% chanel. So all signal reconstruct point in one direction of listening.
% Input var: matPondahf = Ponderation matrix with indice (nbAzimut,nbMicro, nbFreq)
%            matFFT     = matric of spectrogram (nbFreq, nbTimefft)
%            azimut = vecteur of azimut
%            spec =  sepctogram parameter
% Debug:
matPondahf = MAT_POND_vs_azim_h_freq;
matFFT= MAT_ffts_vs_f_h_INT_INT;
%azimut = vec_azimut;
%freq = vec_f_INT;
%save('floData', 'nbV','matPondahf' , 'matFFT','azimut','freq', 'spec')
%load floData.mat

% Default parameter
showFig = [0];


% Indice
indV = ceil(length(azimut)/nbV) * (1:nbV);
azimutV = azimut(indV);

freqTheo = (0:1:spec.Ns-1)*spec.Fs/spec.Ns;
indF = find( (freqTheo >= spec.fmin)&(freqTheo <= spec.fmax));

% Initialize
recon=zeros(spec.Ns,nbV);
for u = 1 : length(indV)
    % Initiate
    matFFTPond = zeros(1,floor(spec.Ns/2)+1);
    % Calcul the fft for the direction and add those freq values to the null matrix
    matFFTPondSum = ( sum(matFFT .* squeeze(matPondahf(u,:,:))',2 ))';
    
    % Add the conj
    matFFTPond(indF) = matFFTPondSum;
    matFFTPond2 = [ matFFTPond conj(fliplr(matFFTPond(2:length(matFFTPond)-1))) ];
    matFFTPondV(:,u) =  matFFTPond2';
    % Reconstruct WAV
    recon(:,u) = ifft(matFFTPond2);
    
end


% calcul des spectrogrammes des voies
for u = 1 : nbV
    [timeV, freqV, matCompV, matFFTV_dB] = COMP_STFT_snapshot(recon(:,u),0, spec.Fs, spec.winSz, spec.rec, spec.wpond, spec.zp);
    
    if u==1 % Get size
        indFr = find( ((freqV>= spec.fmin)&(freqV <= spec.fmax)));
        reconFFT = zeros(nbV,length(timeV),length(indFr));
    end
    reconFFT(u,:,:)=matFFTV_dB(:,indFr);
end

% Return only freq between fmin and fmax
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

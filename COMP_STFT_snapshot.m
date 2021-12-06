function [vec_temps, vec_freq, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB] = COMP_STFT_snapshot(vec_s, t0, fe, Lfft, REC, w_pond, fact_zp)

Nd1_rel=Lfft/2+1; % indice de début
Nf1_rel=length(vec_s)-Lfft/2;% indice de fin
Pas_c=floor(Lfft*(1-REC)); % pas entre deux fenetres FFT pour assurer un taux de recouvrement == REC
vecteur_indice_centre=Nd1_rel:Pas_c:Nf1_rel;
freq_dsp=(0:fact_zp*Lfft-1)*fe/(fact_zp*Lfft);
indi_int = find( (freq_dsp>=0) & (freq_dsp<fe/2));
vec_freq =freq_dsp(indi_int);
vec_temps = vecteur_indice_centre*1/fe + t0;
%h = waitbar(0,'CALCUL DU SPECTROGRAMME IN PROGRESS');
for u = 1 : length(vecteur_indice_centre) % on parcourt les snapshots
    %waitbar(u/length(vecteur_indice_centre))
    ss_loc  = vec_s(vecteur_indice_centre(u)-Lfft/2+1:vecteur_indice_centre(u)+Lfft/2); % on vient sélectionner le snapshot interessant
    a = size(ss_loc);
    aa = size(w_pond);
    if a(1,1)~=aa(1,1)
        w_pond=w_pond';
    end
    fft_loc = fft(ss_loc.*w_pond,fact_zp*Lfft);
    fft_loc_int = fft_loc(indi_int);
    
    if iscolumn(fft_loc_int)
        MAT_t_f_STFT_complexe(u,:)=fft_loc_int';
    else
        MAT_t_f_STFT_complexe(u,:)=fft_loc_int;
    end
end
%close(h)

MAT_t_f_STFT_dB = 10*log10(1/(Lfft*fe)*abs(MAT_t_f_STFT_complexe).^2);
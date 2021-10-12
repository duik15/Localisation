function MAT_POND_vs_azim_h_freq = function_ponderation_azim_h_freq(vec_azimut, xc, yc, vec_f,c)

xc_rel = xc-mean(xc);
yc_rel = yc-mean(yc);

MAT_POND_vs_azim_h_freq = zeros(length(vec_azimut),length(xc_rel),length(vec_f));
for u = 1 : length(vec_azimut)
%     u
%     length(vec_azimut)
    for v = 1 : length(xc_rel)
        for w = 1 :length(vec_f)
    taux =-xc_rel(v)*cos(vec_azimut(u))+-yc_rel(v)*(sin(vec_azimut(u)));
    MAT_POND_vs_azim_h_freq(u,v,w) = exp(sqrt(-1)*2*pi*1/(c/vec_f(w))*taux);
        end
    end
end

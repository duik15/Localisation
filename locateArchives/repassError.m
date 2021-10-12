i=1

size(MAT_s_vs_t_h_INT)
size(MAT_s_vs_t_h)
size(MAT_s)
    
    
ind{i}=  indi_f;
MM1{i} = MAT_s;
MM{i} = MAT_ffts_vs_f_h_INT_INT;
MM2{i}= MAT_ffts_vs_f_h_INT;
MM3{i}  =MAT_s_vs_t_h_INT;
MM4{i}  =MAT_s_vs_t_h;
vs{i} = vec_sfft;
azi{i} = vec_azimut*180/pi;
E{i} = Energie_dB;

%%MAT_ffts_vs_f_h_INT
MM22{i} = nan(size(MM3{i}));
    for u = 1 : Nc
        MM22{i}(:,u)=fft(MM3{i}(:,u));
    end
    

%%
outName = 'pass2'%'OffSet23_voc3_Ns2e14_150-200hz_AvCircle';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%%
i= 2

size(MAT_s_vs_t_h_INT)
size(MAT_s_vs_t_h)
size(MAT_s)

ind{i}=  indi_f;
MM1{i} = MAT_s;
MM{i} = MAT_ffts_vs_f_h_INT_INT;
MM2{i}= MAT_ffts_vs_f_h_INT;
MM3{i}  =MAT_s_vs_t_h_INT;
MM4{i}  =MAT_s_vs_t_h;
vs{i} = vec_sfft;
azi{i} = vec_azimut*180/pi;
E{i} = Energie_dB;

MM22{i} = nan(size(MM3{i}));
    for u = 1 : Nc
        MM22{i}(:,u)=fft(MM3{i}(:,u));
    end

% 
%%

aa = MM1;
figure(6)
plot(aa{1}(:,1))
%xlim([2000, 3000])
xlim([16400, 16900])
hold on
plot(aa{2}(:,1))
xlabel('X')
ylabel('Y')
title('MAT_s_vs_t_h_INT')
legend('Pass1','Pass2')
%M3print(['/Users/Administrator/Documents/MPO/BRing/Data/results/retreive/repass/MAT_s_vs_t_h_INT.png'],'-dpng','-r150')
print(['/Users/Administrator/Documents/MPO/BRing/Data/results/retreive/repass/MAT_ffts_vs_f_h_INT.png'],'-dpng','-r150')

%%
bb = ind;
figure(6)
plot(bb{1})
hold on
plot(bb{2})
xlabel('Azimute')
ylabel('Unkknown')

%%

figure(6)
plot(azi{1},E{1})
hold on
plot(azi{2},E{2})
xlabel('Azimute')
ylabel('Energy')

print(['/Users/Administrator/Documents/MPO/BRing/Data/results/retreive/repass/energy_fct_azimut.png'],'-dpng','-r150')

%%
save('saveData1','MM','MM2','MM22','MM3','azi','E')

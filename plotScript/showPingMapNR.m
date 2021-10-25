%extraire les positions des arrays & pings à partir des .m de Kevin
%6 oct 2021

repwork = 'C:\Moi\Planif2021\Reseaux_circ\Rapport_Technique\Kevin\';
cd (repwork)
lignecote = 'C:\Moi\Cartes\Cote\Golfe_1\COTESUD.DAT'; 

load(lignecote);
%cote = table(lignecote);
% dataID option :  cav = cricle around Anse-a-Valleau
%                  ccd = Cicle around Cloridorme
%                  cmb =   "  " Malbay 
%                  cpc = " " perce
%                  perph = perpendicular transect for Honguedo 
%                  para10h = parallele 10 mn Honguedo
%                  para35h = parallele 35 mn Honguedo
%                  perpp = perpendiculaire Perce
%                  para1p = parallele 1 Perce
%                  para2p = parallele 2 Perce
dlon = 0.02;

trajet = []; tserie = [];

%%% Positions des reseaux

locav = getArrLoc('av');% Anse a Valleau
loccd = getArrLoc('cd');%  Cloridorme
locmb = getArrLoc('mb');%  Mal Bay
locpc = getArrLoc('pc');%  Perce

%loader tous les pts ping
[cAV , tcAV] = getPingLoc('cav'); %cercle Anse a Valleau
[cCL , tcCL] = getPingLoc('ccd'); %cercle Cloridorme
[cMB , tcMB] = getPingLoc('cmb'); %cercle Mal Bay
[cPE , tcPE] = getPingLoc('cpc'); %cercle Mal Bay

[perpH , tperpH] = getPingLoc('perph'); 
[para10H , tpara10H] = getPingLoc('para10h'); 
[para35H , tpara35H] = getPingLoc('para35h'); 
[perpP , tperpP] = getPingLoc('perpp'); 
[para1P , tpara1P] = getPingLoc('para1p'); 
[para2P , tpara2P] = getPingLoc('para2p'); 

%monter les vecteurs temps et trajet
tserie = [tcAV;tcCL;tcMB;tcPE];
tserie = [tserie; tperpH;tpara10H;tpara35H;tperpP;tpara1P;tpara2P];

trajet = [cAV;cCL;cMB;cPE];
trajet = [trajet; perpH;para10H;para35H;perpP;para1P;para2P];
[stserie, iserie] = sort(tserie);
strajet = trajet(iserie,:);

%%%%%
%presenter Honguedo

figure
plot(-cAV(:,2),cAV(:,1),'.b')
grid on; hold on;
plot(-locav(2),locav(1),'or')
plot(COTESUD(:,1),COTESUD(:,2),'k')
%axis([-64.65 -64.45 49.03 49.13])
text(-locav(2)+dlon,locav(1),'Anse-à-Valleau')
plot(-perpH(:,2),perpH(:,1),'*b')
plot(-para10H(:,2),para10H(:,1),'*g')
plot(-para35H(:,2),para35H(:,1),'*k')
%figure
plot(-cCL(:,2),cCL(:,1),'.b')
grid on; hold on;
plot(-loccd(2),loccd(1),'or')
text(-loccd(2)+dlon,loccd(1),'Cloridorme')
%plot(COTESUD(:,1),COTESUD(:,2),'k')
%axis([-65 -64.35 49.03 49.25])
axis([-66.5 -63.5 48.8 50])
title('Honguedo')
print([repwork 'Honguedo'],'-dpng')

%%%%%
%presenter Perce

figure
plot(-cMB(:,2),cMB(:,1),'.b')
grid on; hold on;
plot(-locmb(2),locmb(1),'or')
plot(COTESUD(:,1),COTESUD(:,2),'k')
%axis([-64.65 -64.45 49.03 49.13])
text(-locmb(2)+dlon,locmb(1),'Mal-Bay')
plot(-perpP(:,2),perpP(:,1),'*b')
plot(-para1P(:,2),para1P(:,1),'*g')
plot(-para2P(:,2),para2P(:,1),'*k')
%figure
plot(-cPE(:,2),cPE(:,1),'.b')
grid on; hold on;
plot(-locpc(2),locpc(1),'or')
text(-locpc(2)+dlon,locpc(1),'Percé')
%plot(COTESUD(:,1),COTESUD(:,2),'k')
%axis([-65 -64.35 49.03 49.25])
axis([-64.8 -63.2 47.8 49.2])
title('Percé')
print([repwork 'Perce'],'-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%

%vérifier la chrono
figure
plot(-strajet(:,2),strajet(:,1),'*-')
grid on; hold on;
plot(COTESUD(:,1),COTESUD(:,2),'k')
plot(-strajet(1,2),strajet(1,1),'or')



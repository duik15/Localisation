function [sevenm lstdetect]=wavvisu_multif(fich,dt,tlim,varargin)

% % % % % % % % % % % % % % % % % % %
% ADAPTATION DE ORIGINAL: belugabinvisu_multiff 
% 
% function sevenm=wavvisu_multif(fich,dt,tlim,fichbin,varargin)
%
%
% On s inspire de wavvisu et revisdetect; mais on adapte pour afficher le spectro original et la binarisation par Autotrait_xxxx
%
%
% DESCRIPTION
%        tentative de faire de revisdetect une fct pour visualiser fichier simplement.
%
%
% STRUCTURE ET FONCTIONNEMENT
%
%
% INTRANTS
%           fich:       chemin vers fichier(ch. de caract. ou cellule de ch. de caract.)
%           dt:         largeur de fenetre, (en sec.)
%           tlim:       intervalle de temps a afficher(par sections de dt secondes). mat nfich x 2
%                       si []: montre tout le fichier.
%
% EXTRANTS
%           sevenm: structure qui compile les evenements identifies par l'usager
%                       chaque champs represente une categorie
%                                     et contient une mat n x 3 de forme  [ik, tevenm, fevenm]
%                                               ik:         les indices des spectrogrammes dans lequel se trouve l evenement.
%                                               tevenm:     temps(en sec) de l evenm, a partir du debut du fichier acoustique(et non du debut du segment.)
%                                               fevenm:     freq.(en Hz) de l evenm.
%
%                               champs par defaut: 'voc':       la vocalise recherchee
%                                                  'bruit':     bruit confondant:  offre la possibilite d identif. les bruits
%                                                                                  les plus susceptible de poser probleme
%                                                  'indecis':
%                                                  'saute':     spectrogramme n ayant aucun evenement(vestige revisdetect, possiblement a supprimer.)
%
%                               champs donnes dans champs 'categ' de struct optionnelle.(voir section OPTIONS)
% 
%           lstdetect:       cellule des segments acoustiques (dim N x 2) correspondant a chaque fenetre de visualisation et classification.
%                            1re colonne nom fich wav(chemin complet),
%                            2e colonne, temps de centre du segment, en seconde, a partir du debut du fichier
%
%   OPTIONS:    l ensemble des options sont inclues dans une structure 'param' (passe en varargin)
%               struct 'param' sera passee a wavvisuconstruct.
%
%
%
%           param:       structure contenant les champs qui fixent les parametres de construction des spectrogrammes
%                  champs obligatoires:
%                                        'FSd':             freq. de reechantillonage de signal
%                                        'frame':           Largeur de fen�tre(en nb de points) pour calculer le spectre des fichiers wav.
%                                        'padding':         (zero-padding): augmente le nb de freq du spectro par: nfft=frame*padding;
%                                        'chevauchement':   Recouvrement du fenetrage(entre 0 et 1).
%                                        'fmin':            freq max et min affichees (initialement concue pour correspondre au freq max et min du noyau de correlation spectrale utiliser pour la detection)
%                                        'fmax'
%                                        'dtvoc'            intervalle de temps affiche avant et apres chaque detection(on construit donc le spectrogramme
%                                                           d une duree de 2*dtvoc (arrondie au nb entier de pas de temps pres)
%
%
%                   champs optionnels
%                                         'profil':        selon le profil, on fixe les valeurs par defaut adaptees au profil.
%                                                          pour instant:
%                                                          'rorqual':   pour rorqual bleu et commun
%                                                          'beluga':    pour beluga(par defaut)
%
%                                         'dtcontxt'      intervalle de temps affiche avant et apres chaque detection(elargi pour apercu contexte)
%                                                         (on construit donc le spectrogramme d une duree de 2*dtcontxt (arrondie au nb entier de pas de temps pres)
%
%                                          'audio'        si =0 ou champs absent: pas d audio
%                                                             si =1, pour chaque detect., construit objet audio:    de la vocalise (si sans contexte)
%                                                                                                                 du contexte (si avec contexte)
%
%                                                         retourne 'haud': cellule d etiquettes(handles) vers les objets audio (dim:  1 x N ou N x 1(si sans contexte))
%                                                                                                                              (dim:  N x 2 (si avec contexte) )
%                                                                                                                                   colonne 1: etiquettee vers obj voc.
%                                                                                                                                   colonne2: etiquette vers obj contexte
%
%
%                                                          on pourra utiliser les etiquettes pour jouer les fichiers wav(reechantillones a freq FSd).
%
%                                          'amplivol'    option: amplification du volume audio.
%                                                               si =0 ou champs absent: garde signal original
%                                                               si =1, amplifie le signal: amplif maximalement le signal(le max(signal)=+-1)
%                                                                                   (divise par max(signal); donc signal maximal, sans saturation artificielle))
%
%
%                                           'fsmin'      force freq minimum de echantillonage, selon la carte de son utilisee
%                                                        si FSd<fsmin et option audio, force l audio a cette freq minimum. (spectro reste construit a partir Wav a FSd)
%                                                        (par defaut: ifsmin=1(utilise l option) et fsmin 4000Hz)
%                                                                     un meilleur defaut serait   ifsmin=0; car plus rapide a exectution
%
%
%                                           'flim'        option permettant d etendre la visualistion en frequence ds le contexte(seulement utile si utilise option de contexte)
%                                                               param.flim=[fmincontxt fmaxcontxt]; vecteur avec limites min. et max. de freq desirees.
%                                                                                                     si freq n etendent pas la plage p/r a [fmin fmax] on les ignorent
%
%                                           'typefen'    option permettant
%                                           choisir le type de fenetrage
%                                           pour le spectrogram (defaut 'hanning')
%                                           (voir help window)
%
%                                            'msqglob'    option: affichage du msqglob de specbin
%                                                               si =0 ou champs absent: affiche pas
%                                                               si =1, affiche
% 
%                                            'msqexterne' option: affichage de msqexterne, un msq qui se trouve ds fichier autre que fichbin
%                                                                 il est combine avec msqglob si celui-ci est actif
%                                                          sous-struct de champs:
%                                                           (...).msqexterne.drapeau: drapeau logique; presence ou non de msqexterne
%                                                           (...).msqexterne.fichmsq: chemin complet vers fichier de msq externe; correspondant a chaque fich et fichbin(cellule de ch. de caract.);
%                                                           (...).msqexterne.varmsq:  nom de la variable du msqexterne dans le fichier(ch.caract.)
%   
% 
%                                           'log'    option log: fichier d annotation de revision
%                                                    structure, de ch.:
%                                                    param.log.flog:  nom du fichier log(chemin complet)
% 
%                                           'icanal' canal a charger ds fichwav (defaut: 1)                                          
% 
%                         param.log.dt:    largeur de la fenetre de visualisation (en fait, largeur du pas de temps entre 2 visualisation, car les fenetres se chevauchent un peu.)
% 
%                                          -cree(ou ) un fichier d annotations csv contenant nom du fich vu, avec temps ti et tf, puis une colonne pour les notes
% 
%                                          -prend des notes de 3 manieres: 1- le classement en voc ou indecis est consigne(de meme que les categ.supp.) (si par contre, choisi 'saute', 'refuse' ou une scategsupp, retrouve pas ds log)
%                                                                          2- action 'l': entraine apparition d une fenetre ou ecrire ses notes
%                                                                          3- action 'k': entraine apparition d une fenetre ou tape des raccourcis a une lettre pour des
%                                                                                         annotations frequentes.
%                                                                                          par exemple: les lettres: 'mb' ajoutent ds le log, '--Beluga manquee--Bateau
%                                                                                          binarise'
% 
%                                                                     les codes sont:          {'v','b','c','a','m'};
%                                                                     et correspondent aux annotations   {'Beluga interess','Bateau binarise','Cognement binarise','Autre animal','Beluga manque'};
% 
%
%
% EXEMPLE: appel en modifiant les param par defaut:
% wavvisu('B8A80042.WAV',60,[15,24]*60,struct('accelaud',1,'fmin',1000,'fmax',5000,'FSd',10000,'chevauchement',0.125,'frame',1024))
%
% % % % % % % % % % % % % % % % % % % %
% % HISTORIQUE DE MODIF, SAM
% -2016sept30:  retabli option de contxt qui marchait pas.(manquait le if isfield(param,'dtcontxt'))
% -2017oct20:	ajuste les profils, chmps defaults: frame -> nfen; FSd-> frs.
%         							ajoute defaults.tc=0 pour tous profils.
% -2018dec17:   profil beluga, change defaults.fmin=eps vers defaults.fmin=100
%                              frequence a peine audible, mais eps peut avoir effet de rendre signal nul apres filtration pssbnd.                              
% -2020fev13:   option msqglob: pour affiche pas de temps masques du specbin
% 
% -2020fev17:  option log: pour fichier d annotation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% auteur: Samuel Giard
% date initiale:
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% % auteur:          samuel giard
% % date initiale:    2020nov18


warning('SG 2021jan19: incertain sur état d adaptation de cette fct!!! Possible que soit pas encore complète; a tester!')


% % % % % % % % % % % % % % % % % % %
% % % ajout de parametres de revisio


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % %   !!!!modifier  selon les besoins!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% % % % % % % % % %
% % struct d option

% % %  % decommenter et ajuster pour activer l option

% % %       %option categ. supp.: construit structure de categories supplementaire(pour les touches reservees, voir help revisdetect, section STRUCTURE ET FONCTIONNEMENT(2))
%     option.scategsupp=struct('categ',{'bosse','autre'},'touche',{'b','a'});
%     option.contxt=1;        %option affichage du contexte
option.db=1;              %option spectrogram en dB (1:dB, 0:lineaire)
%!!!!!!!! fin modif !!!!!!!!!!!

 profil='rorqual'; %defaut
%profil='beluga'; %defaut

if ~isempty(varargin)
    
    if  length(varargin)~=1  || ~isstruct(varargin{1})
        error('wavvisu:varargin',['Arguments opionnels invalides: on accepte uniquement une structure qui definit l ensemble des options'])
    else
        
        if isfield(varargin{1},'profil')
            profil=varargin{1}.profil;
        end
    end %if ||
end %varargin


%valeurs par defaut valides pour tous les profils
defaults.tc=0;
defaults.msqglob=false;
defaults.msqexterne.drapeau=false;
defaults.log.flog='';            % pas de fich log
defaults.log.dt=0;            % pas de fich log
switch profil
    case 'rorqual'
        
        %%%%%%   champs obligatoires %%%%%%%%%%%%%%%%%%%
        
        % % % % % % % % % % % % % % % % % %
        % % % definition des champs de structure de parametre
        
        %defaults.FSd=500;             % frequences de reechantillonage
        %defaults.frame=512;            % fenetre de fft(nb de points)
        defaults.FSd=500;             % frequences de reechantillonage
        defaults.frame=512;            % fenetre de fft(nb de points)
        %         defaults.nfft=9*512;             % nb de freq ds fft
        
        defaults.padding=9;            % bourrage de 0 (zero-padding)
        defaults.chevauchement=0.95;  % chevauchement de fenetres
        defaults.fmin=5;              % freq. minimum de l affichage
        defaults.fmax=100;             % freq. maximum de l affichage
        
        defaults.accelaud=10;   %   facteur d acceleration du signal.
        
        % % % % % % % % % % % % % %
        
        % % comme ca, les segment s entrecoupent un peu.(1.9 au lieu de 2)
        defaults.dtvoc=dt/1.9;    %intervalle de temps affiche avant et apres chaque detection(en sec.)
        
        %%%%%% champs optionnels %%%%%%%%%
        
        % defaults.dtcontxt=15;   %intervalle de temps affiche avant et apres chaque detection(elargi pour apercu contexte)(en seconde)
        defaults.audio=1;       % si =0 ou champs absent: pas d audio
        %                      % si =1, pour chaque detect., construit objet audio:    de la vocalise (si sans contexte)
        %                                                                              du contexte (si avec contexte)
        defaults.amplivol=1;      % 'amplivol'     si =0 ou champs absent: garde signal original
        %                    si =1, amplifie le signal: (divise par max(signal); donc signal maximal, sans saturation artificielle))
        
        %
        % defaults.flim=[50 300];  %  option permettant d etendre la visualistion en frequence ds le contexte(seulement utile si utilise option de contexte)
        %                                                               defaults.flim=[fmincontxt fmaxcontxt]; vecteur avec limites min. et max. de freq desirees.
        
    
    case 'rorqual'
        
        %%%%%%   champs obligatoires %%%%%%%%%%%%%%%%%%%
        
        % % % % % % % % % % % % % % % % % %
        % % % definition des champs de structure de parametre
        
        %defaults.FSd=500;             % frequences de reechantillonage
        %defaults.frame=512;            % fenetre de fft(nb de points)
        defaults.FSd=500;             % frequences de reechantillonage
        defaults.frame=512;            % fenetre de fft(nb de points)
        %         defaults.nfft=9*512;             % nb de freq ds fft
        
        defaults.padding=9;            % bourrage de 0 (zero-padding)
        defaults.chevauchement=0.95;  % chevauchement de fenetres
        defaults.fmin=5;              % freq. minimum de l affichage
        defaults.fmax=100;             % freq. maximum de l affichage
        
        defaults.accelaud=10;   %   facteur d acceleration du signal.
        
        % % % % % % % % % % % % % %
        
        % % comme ca, les segment s entrecoupent un peu.(1.9 au lieu de 2)
        defaults.dtvoc=dt/1.9;    %intervalle de temps affiche avant et apres chaque detection(en sec.)
        
        %%%%%% champs optionnels %%%%%%%%%
        
        % defaults.dtcontxt=15;   %intervalle de temps affiche avant et apres chaque detection(elargi pour apercu contexte)(en seconde)
        defaults.audio=1;       % si =0 ou champs absent: pas d audio
        %                      % si =1, pour chaque detect., construit objet audio:    de la vocalise (si sans contexte)
        %                                                                              du contexte (si avec contexte)
        defaults.amplivol=1;      % 'amplivol'     si =0 ou champs absent: garde signal original
        %                    si =1, amplifie le signal: (divise par max(signal); donc signal maximal, sans saturation artificielle))
        
        %
        % defaults.flim=[50 300];  %  option permettant d etendre la visualistion en frequence ds le contexte(seulement utile si utilise option de contexte)
        %                                                               defaults.flim=[fmincontxt fmaxcontxt]; vecteur avec limites min. et max. de freq desirees.
        
        
    case 'beluga'
        
        %%%%%%   champs obligatoires %%%%%%%%%%%%%%%%%%%
        
        % % % % % % % % % % % % % % % % % %
        % % % definition des champs de structure de parametre
        
        %defaults.FSd=2^13;             % frequences de reechantillonage
        %defaults.frame=2^9;            % fenetre de fft(nb de points)
        %defaults.frs=defaults.FSd;             % frequences de reechantillonage
        %defaults.nfen=defaults.frame;            % fenetre de fft(nb de points)
        
        defaults.FSd=2^13;             % frequences de reechantillonage
        defaults.frame=2^9;            % fenetre de fft(nb de points)
        
        
%         defaults.nfft=2^9;             % nb de freq ds fft
        defaults.padding=3;            % bourrage de 0 (zero-padding)
        defaults.chevauchement=0.5;  % chevauchement de fenetres
        %         defaults.fmin=800;              % freq. minimum de l affichage
        %         defaults.fmax=3.5e3;             % freq. maximum de l affichage
%         defaults.fmin=eps;              % freq. minimum de l affichage
        defaults.fmin=100;              % freq. minimum de l affichage
        defaults.fmax=4e3;             % freq. maximum de l affichage
        
        defaults.accelaud=1;   %   facteur d acceleration du signal.
        
        % % % % % % % % % % % % % %
        
        % % comme ca, les segment s entrecoupent un peu.(1.9 au lieu de 2)
        defaults.dtvoc=dt/1.9;    %intervalle de temps affiche avant et apres chaque detection(en sec.)
        
        %%%%%% champs optionnels %%%%%%%%%
        
        % defaults.dtcontxt=15;   %intervalle de temps affiche avant et apres chaque detection(elargi pour apercu contexte)(en seconde)
        defaults.audio=1;       % si =0 ou champs absent: pas d audio
        %                      % si =1, pour chaque detect., construit objet audio:    de la vocalise (si sans contexte)
        %                                                                              du contexte (si avec contexte)
        defaults.amplivol=1;      % 'amplivol'     si =0 ou champs absent: garde signal original
        %                    si =1, amplifie le signal: (divise par max(signal); donc signal maximal, sans saturation artificielle))
        
        %
        % defaults.flim=[50 300];  %  option permettant d etendre la visualistion en frequence ds le contexte(seulement utile si utilise option de contexte)
        %                                                               defaults.flim=[fmincontxt fmaxcontxt]; vecteur avec limites min. et max. de freq desirees.
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SECTION PARAMINIT : Initialisation des parametres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %
%    a partir de struct param
%   si le champs requis ne s'y trouve pas, on prend la valeur par defaut definie

if ~isempty(varargin)
    
    if  length(varargin)~=1  || ~isstruct(varargin{1})
        error('wavvisu:varargin',['Arguments opionnels invalides: on accepte uniquement une structure qui definit l ensemble des options'])
    else
        
        param=varargin{1};
        
    end %if ||
else
    param=struct;
end %varargin

cch=fieldnames(defaults); %trouve les champs necessaire a l analyse (doivent exister ds defaults)
nch=length(cch);

for ich=1:nch
    if ~isfield(param,cch{ich})
        %                 disp([cch{ich} 'pas passe en param.'])
        eval(['param.' cch{ich} '=' 'defaults.' cch{ich} ';']);
    end
end %ich


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % %  prelimiaires de revision

% % format requis.cellule
if ischar(fich)
    fich={fich};
end



% % dfn des ch. de acquisition de donnee: 'fs','nech','dureeacoust' par lecture du 1er fich de serie.

% MODIF MULTIFICH
nfich=length(fich);
fs=zeros(nfich,1);
nech=zeros(nfich,1);
dureeacoust=zeros(nfich,1);

for ifich=1:nfich
    [~,fs(ifich,1)]     =       wavread(fich{ifich},1,'native');
    tempn                           =       wavread(fich{ifich},'size');
    nech(ifich)       =       tempn(1);
    dureeacoust(ifich)=       nech(ifich)/fs(ifich);
end %ifich
clear  tempn


% % % % dfn des fenetres de visu:
% % dt: longueur d une fenetre
% %


if isempty(tlim) %tout le fich
    
    
    nseg=floor((dureeacoust-2*param.tc)/dt);
    tvec=cell(nfich,1);
    lst=cell(sum(nseg),2); %nb tot de segments
    
    for ifich=1:nfich
        tvec{ifich}=dt/2+param.tc:dt:(nseg(ifich)-1/2)*dt+param.tc;
%         ou encore
%         tvec{ifich}=linspace(param.tc+dt/2, dureeacoust-param.tc,nseg(ifich));

        disp(['la fin va etre coupee de ' num2str(rem(dureeacoust(ifich)-2*param.tc,dt)) ' secondes.'])
        % keyboard
        sfich= repmat(fich{ifich},nseg(ifich),1);
        
        lst(( sum( nseg(1:ifich-1))+1 : sum(nseg(1:ifich) ) ),1) =  cellstr(sfich);
        % lst(:,2)=num2cell( pas.*(1:nseg)' );
        lst(( sum( nseg(1:ifich-1))+1 : sum(nseg(1:ifich) ) ),2) =  num2cell(tvec{ifich}' );
        
        warning('wavvisu:tlim_vide','ATTENTION: tlim vide; on prend [o dureeacoust]; a verifier!!!(incidence seulement ds affichage je crois.)')
        % tlim=[0,dureeacoust];
    end
    tlim=[zeros(nfich,1),Inf(nfich,1)]; %rendu ici, utilise seulement tlim(1); sinon, on pourra adapter mieux.
    
    lstdetectcompl=lst;
    ndetecttot=size(lstdetectcompl,1);
    
else %partie de fich limit. par tlim.


%     nseg=floor((dureeacoust-2*param.tc)/dt);
    nseg=ceil((diff(tlim,[],2))/dt);
    
    tvec=cell(nfich,1);
    lst=cell(sum(nseg),2); %nb tot de segments
    
    for ifich=1:nfich
%         tvec{ifich}=tlim(ifich,1)+dt/2:dt:(nseg(ifich)-1/2)*dt+param.tc;
%         tvec{ifich}=(floor(tlim(ifich,1)/dt)+1/2)*dt:dt:dt*(floor(tlim(ifich,2)/dt)+1/2);

tvec{ifich}=(tlim(ifich,1) + dt/2 : dt : tlim(ifich,1) + dt*(nseg(ifich)-1/2));        
%         TEMP!!!!
        if length(tvec{ifich})~=nseg(ifich)
            warning('nb tvec ~= nseg; a investiguer.')
            keyboard
            
        end
%         ou encore
%         tvec{ifich}=linspace(param.tc+dt/2, dureeacoust-param.tc,nseg(ifich));

        disp(['la fin va etre coupee de ' num2str(rem(diff(tlim(ifich,:),[],2),dt)) ' secondes.'])
        % keyboard
        sfich= repmat(fich{ifich},nseg(ifich),1);
        
        lst(( sum( nseg(1:ifich-1))+1 : sum(nseg(1:ifich) ) ),1) =  cellstr(sfich);
        % lst(:,2)=num2cell( pas.*(1:nseg)' );
        lst(( sum( nseg(1:ifich-1))+1 : sum(nseg(1:ifich) ) ),2) =  num2cell(tvec{ifich}' );
        
        % tlim=[0,dureeacoust];
%         tlim=[zeros(sum( nseg),1),Inf(sum( nseg),1)]; %rendu ici, utilise seulement tlim(1); sinon, on pourra adapter mieux.
    end
    
    lstdetectcompl=lst;
    ndetecttot=size(lstdetectcompl,1);
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %construction des objets preliminaires a revision

lstdetect=lstdetectcompl;

% cspec=cell(size(lstdetect,1),2);
% cTspec=cell(size(lstdetect,1),2);
% cFspec=cell(1,2);

if param.msqglob
%     warning('Necessaire; a revoir!!!--option msqexterne')
%     if param.msqexterne.drapeau
%         error('On peut appliquer seulement un masque; choisir msqglob ou msqexterne')
%     else
        cmsqglob=cell(size(lstdetect,1),1);
%     end
end

% [cspec(:,1),cTspec(:,1),cFspec(1),haud]=wavvisuconstruct(lstdetect,param);
[cspec,cTspec,cFspec,haud]=wavvisuconstruct(lstdetect,param);
% 
% for ifich=1:nfich
% 
%     load(fichbin{ifich},'Tspecred','Fspec','specbinred')
%     if param.msqglob
%         load(fichbin{ifich},'msqglobred')
%     end
%     
%       
%     if param.msqexterne.drapeau
%         load(param.msqexterne.fichmsq{ifich},param.msqexterne.varmsq)
%         eval(['msqexterne=' param.msqexterne.varmsq ';']);
%         clear(param.msqexterne.varmsq)
%         
%         param.msqglob=true;
%         if exist('msqglobred','var')
%             msqglobred=msqglobred & msqexterne; % combine les msqglob et msqexterne
%         else
%             msqglobred=msqexterne;              % msqexterne seulement, car msqglob existe pas
%         end
%     end
%     
%     
% %     Tspecred=Tspecred + param.tc + tlim(1); % redefini le temps pour annuler la translation imposee par Autotrait.
% Tspecred=Tspecred + param.tc + tlim(ifich,1); % redefini le temps pour annuler la translation imposee par Autotrait.
%     
%     % ajoute le temps coupe. Et la tlim(qui vient de trelsa)
%     %     for ifich=1:
%     if ifich==1
%         cFspec{2}=Fspec;
%     end
%     %
%     % for ik=1:size(cTspec,1)
%     %     it=Tspecred>=cTspec{ik,1}(1) & Tspecred<=cTspec{ik,1}(end); % select. pas inclus ds segment ik de revision
%     %     cTspec{ik,2}=Tspecred(it);
%     %     %        cspec{ik,2}=S50msq2red(:,it);
%     %     cspec{ik,2}=specbinred(:,it);
%     % end
%     
%     for ik=1:nseg(ifich)
%         it=Tspecred>=cTspec{sum( nseg(1:ifich-1))+ik,1}(1) & Tspecred<=cTspec{sum( nseg(1:ifich-1))+ik,1}(end); % select. pas inclus ds segment ik de revision
%         cTspec{sum( nseg(1:ifich-1))+ik,2}=Tspecred(it);
%         if param.msqglob
%             cmsqglob{sum( nseg(1:ifich-1))+ik,1}=msqglobred(it);
%         end
%         cspec{sum( nseg(1:ifich-1))+ik,2}=specbinred(:,it);
%     end
%     
%     % ( sum( nseg(1:ifich-1))+1 : sum(nseg(1:ifich) ) )
% end %ifich

% % % fin construction
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % %  revision manuelle
%modif pour eviter cas ou genere pas audio
if isfield(param,'audio')
    if param.audio==1;
        option.audio=haud;     %option audio: 'etiq': etiquette objet audio
    end
end
%si accelaud, passe facteur a classmanuelle
if isfield(param,'accelaud')
    option.accelaud=param.accelaud;     %option accelaud: facteur d acceleration du signal.
end

if isfield(param,'dtcontxt')
    option.contxt=true;     %option contxt: affich contxt
end

if isfield(param,'tabs')
    option.tabs=param.tabs+sec2jour(tvec);     %option contxt: affich contxt
end

if isfield(param,'scategsupp')
    option.scategsupp=param.scategsupp;     %option categ supplemtr.
end

% if isfield(param,'msqglob')
%     option.msqglob=cmsqglob;     %option msqglob: passe les masques de bruit large bande
% end
% 
% if isfield(param,'log')
%     if ~isempty(param.log)
%         option.log=param.log;
%         
% %         option.log.replog=param.log.replog
% %         option.log.id=param.log.id
% %         option.log.flog=param.log;     % option log: fichier d annotation de revision
%         option.log.dt=dt;         % largeur de la fenetre de visualisation (en faite, largeur du pas de temps entre 2 visualisation, car les fenetres se chevauchent un peu.)
%     end
% end

if isfield(param,'baba')
    option.baba=param.baba; %option affichage limites de bandes baba
end

% error('RENDU A CHANGER PALETTE DE VISUALISATION, COMME DS REVISBELUGA.')

if exist('option','var')
    sevenm=classmanuelle(lstdetect,cspec,cTspec,cFspec,option); %revision manuelle
else
    sevenm=classmanuelle(lstdetect,cspec,cTspec,cFspec); %revision manuelle
end

% % % % % fin revision manuelle
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % REMARQUE POUR L UTILISATION DU RESULTAT DE REVISION
%
%   UN EXEMPLE DE MANIERE DE CREER UNE NOUVELLE VARIABLE CONTENANT QUE LES DETECTIONS CLASSEES DANS UNE CATEGORIE PARTICULIERE:
%
%
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     % % % % compile revision
%
%     %ici, pourrait supprimer cet etape, cause dedoublement de info.
%     % mais pratique
%     % autrement, le recreer chaque fois.
%
%     chrev=fieldnames(srevis);    %liste des champs de option
%     nchrev=length(chrev);
%
%     for ch=chrev'
%         lstrevis.(ch{1})=lstdetect(srevis.(ch{1}), : );
%         cspecrevis.(ch{1})=cspec(srevis.(ch{1}), : );
%         cTspecrevis.(ch{1})=cTspec(srevis.(ch{1}), : );
%         cFspecrevis.(ch{1})=cFspec( : );
%     end %ich
%
%     % % % % fin compile
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %








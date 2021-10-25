function [ploc , time] = getPingLoc(dataID,varargin)
% [ploc, ptime ]  = getPingLoc('cav')
% This fonction return the location and time of any set of ping emission
% For circular ping, the ping number 1 is located the most north
% Times are in UTC.
% dataID option :  cav = cricle around Anse-a-Valleau
%                  ccd = Cicle around Cloridorme
%                  cmb =   "  " Malbay / cpc = " " perce
%                  perph = perpendicular transect for Honguedo 
%   last update 14/09/2021

% Update 5 oct 2021: N.Roy pour ajouter les pings rectilignes
%                  para10h = parallele 10 mn Honguedo
%                  para35h = parallele 35 mn Honguedo
%                  perpp = perpendiculaire Perce
%                  para1p = parallele 1 Perce
%                  para2p = parallele 2 Perce



switch lower(dataID)
    case 'cav'   % -------- Circle at Anse-a-Valleau -----------
        ploc = [ 49.0953 64.5370;
            49.0945 64.5328;
            49.0922 64.5303;
            49.0893 64.5304;
            49.0871 64.5329;
            49.0863 64.5370;
            49.0872 64.5410;
            49.0894 64.5435;
            49.0922 64.5434;
            49.0945 64.5409];
        time=[ datetime(2021,07,15,18,11,00);
            datetime(2021,07,15,18,14,00);
            datetime(2021,07,15,18,16,20);
            datetime(2021,07,15,18,18,58);
            datetime(2021,07,15,18,21,37);  % Was written 18:31:37 in logbook
            datetime(2021,07,15,18,25,15);
            datetime(2021,07,15,18,27,56);
            datetime(2021,07,15,18,30,12);
            datetime(2021,07,15,18,05,02);
            datetime(2021,07,15,18,08,44);];
    case 'ccd'   % -------- Circle at Cloridorme -----------
        ploc = [49.1968 64.8351;
            49.1977 64.8315;
            49.1967 64.8274;
            49.1945 64.8248;
            49.1916 64.8247;
            49.1887 64.8313;
            49.1895 64.8356;
            49.19175 64.838;
            49.1895 64.8273;
            49.1945 64.8378;];
            % 49.1931 64.8314                   % Center array emission
        time=[ datetime(2021,07,15,14,25,50);      % Two emission at same locaion
            datetime(2021,07,15,14,29,00);
            datetime(2021,07,15,14,32,00);
            datetime(2021,07,15,14,34,22);
            datetime(2021,07,15,14,37,00);
            datetime(2021,07,15,14,42,00);
            datetime(2021,07,15,14,46,00);
            datetime(2021,07,15,14,19,17);
            datetime(2021,07,15,14,39,42);
            datetime(2021,07,15,14,22,36) ];
            %datetime(2021,07,15,14,51,00);     % Center array emission         
           
        case 'cmb'   % -------- Circle at Malbay -----------
        % Malbay ping have been made by an knock on the boat
        ploc = [48.6084 64.1863;
                48.6076 64.1824;
                48.6053 64.1800;
                48.6025 64.1799;
                48.6002 64.1823;
                48.5994 64.1863;
                48.6003 64.1903;
                48.6026 64.1928;
                48.6054 64.1926;
                48.6076 64.1903;]; % center emission
                %48.6039 64.1862;];
        time=[ datetime(2021,08,01,14,30,15);
            datetime(2021,08,01,14,32,14);
            datetime(2021,08,01,14,34,23);
            datetime(2021,08,01,14,36,24);
            datetime(2021,08,01,14,38,12);
            datetime(2021,08,01,14,39,55);
            datetime(2021,08,01,14,41,29);
            datetime(2021,08,01,14,43,42);
            datetime(2021,08,01,14,45,00);
            datetime(2021,08,01,14,48,48);];
            %datetime(2021,08,01,14,49,05);  % center emission
        case 'cpc'   % -------- Circle at Perce -----------
        ploc = [48.5356 64.1984;
                48.5348 64.1945;
                48.5326 64.1921;
                48.5298 64.1919;
                48.5275 64.1945;
                48.5267 64.1986;
                48.5277 64.2024;
                48.5298 64.2049;
                48.5325 64.2048;
                48.5349 64.2024;];
                %48.5473 64.1991];
        time=[ datetime(2021,08,01,13,40,37);      % Two emission at same locaion
               datetime(2021,08,01,13,43,17);
               datetime(2021,08,01,13,45,55);
               datetime(2021,08,01,13,47,41);
               datetime(2021,08,01,13,49,34);
               datetime(2021,08,01,13,51,11);
               datetime(2021,08,01,13,53,06);
               datetime(2021,08,01,13,54,44);
               datetime(2021,08,01,13,56,17);
               datetime(2021,08,01,13,57,58);];
               %datetime(2021,08,01,14,00,45);
               
        case 'perph'   % -------- perp Honguedo -----------
        
        ploc = [ 49.1404	64.6809;
            49.1462	64.6811;
            49.1497	64.6776;
            49.1738	64.6566;
            49.2217	64.6179;
            49.3449	64.5195;
            49.4914	64.4122;
            49.6195	64.30];
        time=[ datetime(2021,07,13,23,50,00);
            datetime(2021,07,14,00,25,00);
            datetime(2021,07,14,00,38,42);
            datetime(2021,07,14,01,11,00);
            datetime(2021,07,14,01,58,00);
            datetime(2021,07,14,03,26,59);
            datetime(2021,07,14,05,01,47);
            datetime(2021,07,14,06,21,16)];

        
         case 'para10h'   % -------- parallele 10 Honguedo -----------
        ploc = [ 49.474	65.3528;
            49.3877	65.106;
            49.3258	64.9243;
            49.2836	64.8001;
            49.255	64.7171;
            49.2212	64.6201;
            49.18763	64.52298;
            49.15828	64.43797;
            49.11612	64.31700;
            49.05248	64.13545];

        time=[ datetime( 2021,07,14,20,30,16);
            datetime(2021,07,14,21,47,13);
            datetime(2021,07,14,22,47,00);
            datetime(2021,07,14,23,29,00);
            datetime(2021,07,15,15,08,51);
            datetime(2021,07,15,15,22,21);
            datetime(2021,07,15,15,51,00);
            datetime(2021,07,15,15,47,38);
            datetime(2021,07,15,16,02,20);
            datetime(2021,07,15,16,26,30)];

         case 'para35h'   % -------- parallele 35 Honguedo -----------  
        ploc = [49.6703	65.1916;
            49.5902	64.9448;
            49.5279	64.7596;
            49.4849	64.6413;
            49.4552	64.5581;
            49.3876	64.3597;
            49.3599	64.2744;
            49.3171	64.1536;
            49.2542	63.968];
        time=[ datetime(  2021,07,14,18,27,00);
            datetime(2021,07,14,17,14,00);
            datetime(2021,07,14,16,15,00);
            datetime(2021,07,14,15,33,00);
            datetime(2021,07,14,15,05,00);
            datetime(2021,07,14,14,02,00);
            datetime(2021,07,14,13,27,00);
            datetime(2021,07,14,12,46,00);
            datetime(2021,07,14,11,44,00)];

                   
        case 'perpp'   % -------- perp Percé -----------
       ploc = [  48.56635	64.15064;
            48.56125	64.12397;
            48.56136	64.13170;
            48.55950	64.11190;
            48.48433	63.71638;
            48.44532	63.39751];
      time=[ datetime(  2021,07,30,14,08,00);
            datetime(2021,07,30,14,21,48);
            datetime(2021,07,30,14,31,58);
            datetime(2021,07,30,14,40,20);
            datetime(2021,07,30,18,31,00);
            datetime(2021,07,30,21,21,27)];
            
        case 'para1p'   % -------- para 1 Percé -----------
       ploc = [  49.06618	63.82037;
            48.89348	63.89817;
            48.76353	63.95297;
            48.67728	63.98993;
            48.61655	64.01623;
            48.47808	64.07422;
            48.41713	64.09987;
            48.33085	64.13575;
            48.20155	64.19178;
            48.02733	64.26250 ];       
        time=[ datetime(   2021,07,31,13,58,40);
            datetime(2021,07,31,15,14,13);
            datetime(2021,07,31,16,15,47);
            datetime(2021,07,31,17,02,58);
            datetime(2021,07,31,17,36,58);
            datetime(2021,07,31,18,43,00);
            datetime(2021,07,31,19,15,20);
            datetime(2021,07,31,20,00,26);
            datetime(2021,07,31,21,00,53);
            datetime(2021,07,31,22,18,43)];
        case 'para2p'   % -------- para 2 Percé -----------
       ploc = [  47.96468	63.93440;
            48.13873	63.86063;
            48.26823	63.80698;
            48.35485	63.77113;
            48.41587	63.74563;
            48.55375	63.68690;
            48.61153	63.66237;
            48.69875	63.62435;
            48.82963	63.56987;
            49.00298	63.49382];
       time=[ datetime( 2021,07,31,02,52,58);
            datetime(2021,07,31,04,30,09);
            datetime(2021,07,31,05,46,34);
            datetime(2021,07,31,06,38,48);
            datetime(2021,07,31,07,20,10);
            datetime(2021,07,30,19,23,20);
            datetime(2021,07,31,08,56,26);
            datetime(2021,07,31,09,44,19);
            datetime(2021,07,31,10,45,10);
            datetime(2021,07,31,12,08,11)];

        

    otherwise
        error(['Can''t location ID. Refer to heater.' ])
end

% making sure the lon is negative
function [matPondahf] = makeMatPond(arrID,azimut,Ns,Fs,varargin)
%%
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'fmin'
                fmin = varargin{2};
            case 'fmax'
                fmax = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end

%%

% Get array info
[aPos arr]  = getArrInfo(arrID);

% Get array needed
azimutR= pi* azimut/180;  % Vector azimut
freq = (0:1:Ns-1)*Fs/Ns;  % Vector frequence 

% Get freq lim
if ~exist('fmin')
    fmin = min(freq);
end 
if ~exist('fmax')
    fmax = max(freq);
end 

% Create an azimut vecteur  
indF = find( (freq >= fmin)&(freq <= fmax));
freqF = freq(indF);


% Matrice of ponderation
matPondahf = zeros(length(azimutR),length(arr.xh),length(freqF));
Ncapt = length(arr.xh);
for u = 1 : length(azimutR)
    for v = 1 : length(arr.xh)
        d_marche =-arr.xh(v)*sin(azimutR(u))+-arr.yh(v)*(cos(azimutR(u)));
        %MAT_POND_vs_azim_h_freq(u,v,:) = 1/Ncapt*exp(-sqrt(-1)*2*pi*vec_f_INT*d_marche/c);
        matPondahf(u,v,:) = 1/Ncapt*exp(-sqrt(-1)*2*pi*freqF*d_marche/arr.c);
    end
end

end % Function

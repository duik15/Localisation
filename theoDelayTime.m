%  Create theoritacal time delay for circular array
clear all

% Default parameter
sv =  1490; % Sound Speed


% Configuration values
sLoc = [49.1917 64.8380];

% Hydrophone location
arrLoc = getArrLoc('cd');
hLoc = createHLoc(arrLoc);
nbH = length(hLoc(:,1));

% Compute distance beetween source and hydrophone
for ii=1:nbH    
[dSH(ii), aSH(ii), aHS(ii) ] = m_idist(sLoc(2),sLoc(2),hLoc(ii,2),hLoc(ii,1));
end


% Pair indice
nH = [1:10];
n = 1;				% time-delay number
for ii = 1 : nbH
  for jj = ii+1 : nbH
    m1(n) = nH(ii);
    m2(n) = nH(jj);
    n = n + 1;
  end
end


% Time to reach the hydrophone
tSH  = dSH / sv;

% Delay time
for ii = 1: length(m1)
   dt12(ii) = tSH(m2(ii)) - tSH(m1(ii)) ;   
end



%% Figure 1

figure(1)
plot([2:nbH],dt12(m1==1))
ylabel('Delay (s)')
xlabel('Hydrophone')

print('-dpng','-r150','results/delayTimeTho.png')




function  machinName = fgetMachin()
% This script verfy is we are working on a Local machin or mingan
pwdname = pwd;

    if strcmp(pwdname(1:6), '/share')
        machinName ='mingan';
    elseif strcmp(pwdname(1:6), '/Users')
        machinName ='mac';
    else
        error ('Wrong path ! Probably in home! ')
    end
    
end
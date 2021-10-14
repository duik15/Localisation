function  machineName = getMachine(varargin)
% This is use to get information about on which computer the user is
% working and send information of path related to this machine
% Add your as needed

pwdname = pwd;

    if strcmp(pwdname(1:6), '/share')
        machineName ='mingan';
    elseif strcmp(pwdname(1:6), '/Users')
        machineName ='mac';
    else
        error ('Wrong path ! Add you machine to getMachine! ')
    end
    
end
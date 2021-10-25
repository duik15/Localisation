function  machineName = getMachine(varargin)
% This is use to get information about on which computer the user is
% working and send information of path related to this machine
% Add your as needed

pwdname = pwd;
    


    if strcmp(pwdname(1:6), 'C:\Use')
        % Here are windows computer
        
        % Get computer name
        [~, compName]  = system('hostname');
        
        if length(compName > 14) && strcmp(compName(1:14), 'WLQCIML0115503')
            machineName ='dellKev';
        else
            % Add computer 25 or your own computer
            error ('Wrong path ! Add you machine to getMachine! ')
        end
    elseif strcmp(pwdname(1:6), '/Users')
        % Here are mac computer
        machineName ='mac';
    else
        error ('Wrong path ! Add you machine to getMachine! ')
    end
    
end
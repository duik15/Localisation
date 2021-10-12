% This managed the input argument of locateSound

while ~isempty(varargin)
        switch lower(varargin{1})
            case 'outname'
                outName = varargin{2};
            case 'folderin'
                folderin = varargin{2};
            case 'folderout'
                folderout = varargin{2};
            case 'showfig'
                showFig = varargin{2};
            case 'savedata' 
                saveData = varargin{2};
            case 'printfig'
                printFig = varargin{2};
            case 'nbPk' 
                nbPk = varargin{2};
            case 'ns'
                Ns = varargin{2};
            case 'buffer'
                buffer = varargin{2};
            case 'fmin' 
                fmin_int = varargin{2};
            case 'fmax'
                fmax_int = varargin{2};
            case 'moreInfo'
                moreInfo = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end
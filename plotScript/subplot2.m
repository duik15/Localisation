function [ positions, meshgrid2] = subplot2(varargin)
%Return postion of subplot with more option tahtn subplot
%sub_pos=subplot2(nbx,nby,varargin);
%sub_pos=subplot2(sp,varargin); Where sp is a structure with the option
%
% List of options: height, width, spacex, spacey, fracy = [0.2 0.8],
% fracx,ledge,redge,bedge,tedge


%% Varagin

meshOrder = 'y';
while ~isempty(varargin)
    if isstruct(varargin{1})
        sp = varargin{1};
        varargin(1)=[];
    elseif isnumeric(varargin{1})
        sp.nbx =  varargin{1};
        sp.nby =  varargin{2};
        varargin(1:2)=[];
    else
        switch lower(varargin{1})
            case 'height'
                sp.height=varargin{2};
            case 'width'
                sp.width= varargin{2};
            case 'spacex'
                sp.spacex=  varargin{2};
            case 'spacey'
                sp.spacey=  varargin{2};
            case 'fracy'
                sp.fracy=  varargin{2};
            case 'fracx'
                sp.fracx=  varargin{2};
            case 'meshorder'
                meshOrder = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
        varargin(1:2)=[];
    end

%%
% Set default value if not specified

if ~exist('sp')
    sp = struct();
end
if ~isfield(sp,'height')
    sp.height=20;
end
if ~isfield(sp,'width')
    sp.width=30;
end
if ~isfield(sp,'ledge')
    % we assume that if none edge is specified we need all of them
    sp.ledge=2; sp.redge=3; sp.tedge=1.2; sp.bedge=2;
end
if ~isfield(sp,'spacex')
    sp.spacex=0.5;
end
if ~isfield(sp,'spacey')
    sp.spacey=0.5;
end
if ~isfield(sp,'fracy')
    sp.fracy = repmat(1/sp.nby,1,sp.nby);
end
if ~isfield(sp,'fracx')
    sp.fracx = repmat(1/sp.nbx,1,sp.nbx);
end


% Calcul the size of each panel
subxsize=((sp.width-sp.ledge-sp.redge-sp.spacex*(sp.nbx-1.0)) ) * sp.fracx;
subysize=((sp.height-sp.tedge-sp.bedge-sp.spacey*(sp.nby-1.0)) ) * sp.fracy;

% Calcul the posotion
for i=1:sp.nbx
    for j=1:sp.nby
        
        if i==1
            xfirst= sp.ledge;
        else
            xfirst= sp.ledge+(i-1.0)*sp.spacex +  sum(subxsize(1:i-1));
        end
        
        if j==sp.nby
            yfirst=sp.bedge; % nby-j insted of j-1
        else
            yfirst=sp.bedge+(sp.nby-j)*sp.spacey + sum(subysize(j+1:sp.nby)); % nby-j insted of j-1
        end
        
        positions{i,j}=[xfirst/sp.width yfirst/sp.height subxsize(i)/sp.width subysize(j)/sp.height];
        
    end
end


% Meshgrid for incrematation
if meshOrder == 'y'
meshgrid2 = reshape(1:sp.nbx*sp.nby,sp.nbx,sp.nby)';
else
    meshgrid2 = reshape(1:sp.nbx*sp.nby,sp.nby,sp.nbx);
end

%{
% Save
subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;
    for i=1:nbx
       for j=1:nby
           if j == 1
           yfirst=bottommargin+(nby-j)*(subysize+spacey)+spaceyOne; % nby-j insted of j-1
           else
           yfirst=bottommargin+(nby-j)*(subysize+(spacey)); % nby-j insted of j-1
           end
           xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
           
           if j==1
               sub_pos{i,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth (2*subysize)/plotheight];
           else
                sub_pos{i,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
           end
          
       end
    end
%}


end

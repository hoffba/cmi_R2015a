% CMIclass function
function [MF,labels] = calcMF(self,varargin)
% Calculate Minkowski Functionals
% Either called by CMI GUI or through command line input:
% Syntax: cmiObj.calcMF(MFdata,MFval, ... [optional]'Name'/Value pairs)
%                                    'mode', any of {'<','>','>=','<=','=='}
%                                    'win', [1x3] vector of integers
%                                    'voi', true/false
%                                    'defval', true/false

MF = []; labels = {};

% Initialize defaults
MFdat = self.vec;
MFval = [];
MFmode = '>';
MFwin = inf(1,3);
MFvoi = true;
MFdefval = false;

% If a GUI call, ask for user input:
if (nargin==1) || ((nargin==3) && isa(varargin{1},'matlab.ui.container.Menu'))
    str = {'MF Data ("VOI", "PRM", or image index)', num2str(MFdat);
           'Values (for binarization)',              num2str(MFval);...
           'Mode (p if percentiles, then <,>,>=,<=,==)',MFmode;...
           'Window Size ("inf" for global analysis)',  num2str(MFwin);...
           'Apply VOI to analysis? (y/n)',              'y';...
           'Non-VOI value (1/0)',                       num2str(MFdefval)};
    answer = inputdlg(str(:,1),'Minkowski Functionals',1,str(:,2));
    if isempty(answer)
        return;
    else
        MFdat = str2double(answer{1});
        if isnan(MFdat)
            MFdat = answer{1};
        end
        MFval = str2num(answer{2});
        MFmode = answer{3};
        MFwin = sscanf(answer{4},'%f');
        MFvoi = strcmpi(answer{5},'y');
        MFdefval = logical(str2double(answer{6}));
    end
else
    p = inputParser;
    addRequired(p,'data',@(x)(ischar(x)&&any(strcmpi(x,{'voi','prm'})))...
                           ||(isnumeric(x)&&(x>0)&&(x<self.img.dims(4))));
    addRequired(p,'val',@isnumeric)
    addParameter(p,'mode',MFmode,@(x)ischar(x)&&ismember(x,{'<','>','>=','<=','=='}));
    addParameter(p,'win',MFwin,@(x)isvector(x)&&(length(x)==3));
    addParameter(p,'voi',MFvoi,@islogical);
    addParameter(p,'defval',MFdefval,@islogical);
    parse(p,varargin{:});
    
end

% Minkowski functional analysis:
[MF,labels] = self.img.calcMF(MFdat,MFval,MFwin,...
                              'ApplyMask',MFvoi,...
                              'tmode',MFmode,...
                              'defVal',MFdefval);
                          
% Return results to base workspace and display:
dout = array2table(MF,'VariableNames',labels,...
    'RowNames',cellfun(@num2str,num2cell(MFval),'UniformOutput',false));
disp(dout)
assignin('base','MFresults',dout);


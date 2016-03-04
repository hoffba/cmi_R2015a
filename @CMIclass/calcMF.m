% CMIclass function
function [MF,p] = calcMF(self,varargin)
% Calculate Minkowski Functionals
% Either called by CMI GUI or through command line input:
% Inputs: Name/Value pairs for minkowskiFun.m

% Determine inputs:
if nargin>1
    
end

% Minkowski functional analysis:
[MF,p] = self.img.calcMF(self.vec);
                          
if ~isempty(MF) && ismatrix(MF)
    % Return results to base workspace and display:
    dout = array2table(MF,'VariableNames',p.labels,...
        'RowNames',cellfun(@num2str,num2cell(p.thresh),'UniformOutput',false));
    disp(dout)
    assignin('base','MFresults',dout);
end


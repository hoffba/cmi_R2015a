function R = addToReport(R,ch_str,varargin)

import mlreportgen.report.*;
import mlreportgen.dom.*;

if isempty(R) || ~isvalid(R)
    % Initialize the Report with subject ID
    if nargin==2 && (ischar(ch_str) || isstring(ch_str))
        R = Report(['PipelineOutput_',ch_str],'pdf');
        tR = TitlePage;
        tR.Title = 'Pipeline Results';
        tR.Author = ch_str;
        append(R,tR);
        append(R,TableOfContents);
    end
else
    % Add results to Report chapter
    
    switch ch_str
        case 'seg'
        case 'airways'
        case 'vessels'
        case 'unreg'
        case 'reg'
        case 'prm'
        case 'tprm'
    end
end

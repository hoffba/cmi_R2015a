function cmiObj = cmi
go = true;
onum = 0;
curobjs = evalin('base','who(''cmiObj*'')');
if ~isempty(curobjs)
    no = length(curobjs);
    for i = 1:length(curobjs)
        % Clear any that are not valid
        if ~evalin('base',['isvalid(' curobjs{i} ')'])
            evalin('base',['clear ' curobjs{i}]);
            no = no - 1;
        elseif onum==str2double(curobjs{i}(7:end))
            % Make sure instance number is not reapeated
            onum = onum+1;
        end
    end
    if length(curobjs)==no
        onum = num2str(no); % counting starts at 0, so no "+1" necessary
    end
    if no % At least one cmiObj is already open
        answer = questdlg('CMIclass object already exists, create another?');
        if ~strcmp(answer,'Yes')
            go = false;
        end
    end
end
if go
    oname = ['cmiObj' num2str(onum)];
    cmiObj = CMIclass;
    assignin('base',oname,cmiObj);
    set(cmiObj.h.mainFig,'Name',oname)
end


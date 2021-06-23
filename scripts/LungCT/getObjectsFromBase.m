function out = getObjectsFromBase(classname,n)
% Search for objects of class <classname> in the base workspace
% Optional input <n> = max number of objects returned
if nargin<2
    n = 0;
end
baseVariables = evalin('base','whos');
out = cell(0);
for i = 1:length(baseVariables)
    if (strcmpi(baseVariables(i).class, classname))
        outn = length(out)+1;
        out{outn} = evalin('base',baseVariables(i).name);
        if outn==n
            break
        end
    end
end
if n==1 && ~isempty(out)
    out = out{1};
end
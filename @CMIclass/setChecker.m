% CMIclass function
function setChecker(self,hObject,~)
% Display overlay as a checkerboard

if nargin==3 % Called by UIcontrol
    if isempty(self.checkerM)
        % Ask user for checkerboard size
        answer = inputdlg('Checkerboard Size (integer):','Checkerboard',1,...
                            {num2str(self.checkerd)});
        if ~isempty(answer)
            answer = round(str2double(answer{1}));
            if ~isnan(answer) && (answer>0)
                self.checkerd = answer;
                gchk = true;
                set(hObject,'Checked','on');
            end
        end
    else
        gchk = false;
        self.checkerM = [];
        set(hObject,'Checked','off');
    end
else
    gchk = ~isempty(self.checkerM);
end

if gchk % Create checker mask:
    
    % Grab desired size from the current plot:
    d = size(self.hiover.CData);
    
    % Number of tiles:
    td = ceil(d/(2*self.checkerd));
    
    % Generate the checkerboard mask:
    black = false(self.checkerd);
    white = true(self.checkerd);
    tile = [black white; white black];
    I = repmat(tile,td(1),td(2));
    
    % Set the CMIclass property:
    self.checkerM = I(1:d(1),1:d(2));
    
end

if nargin==3
    self.dispUDslice;
end
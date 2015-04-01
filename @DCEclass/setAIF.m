% DCEclass function
function setAIF(self,val)
    % Inputs: none -- use GUI to select AIF
    %         val = number [1...4]
    %               string file name of empirical AIF vector
    %               vector of actual AIF values
    
    ok = true;
    if nargin==1 % Use GUI to select AIF
        opts = {'Load AIF from .mat',...
                'Mouse - 1',...
                'Mouse - 2',...
                'Rat',...
                'Human'};
        [val,ok] = listdlg('ListString',opts,...
                           'SelectionMode','single',...
                           'Name','Select AIF:',...
                           'InitialValue',self.aifType+1);
        val = val-1;
    end
    if ok
        if (numel(val)==1) && (val==0) % use GUI to find file
            [path,fname] = uigetfile('*.mat',...
                                'Select AIF .mat file:',...
                                'MultiSelect','off');
            if path==0
                val = '';
            else
                val = fullfile(path,fname);
            end
        end
        if ischar(val) && ~isempty(val) % load from file
            if exist(fname,'file')
                x = load(fname);
                if isfield(x,'aif')
                    val = x.aif;
                else
                    warning('Invalid AIF file!');
                end
            end
        end
        
        if isnumeric(val) && numel(val)==numel(self.xdata)
            % Assume input values are SI, so need to convert to [Gd]
            self.aif = self.R2gad(self.S2R(val));
            self.aifType = 0;
        elseif numel(val)==1 && ismember(val,1:4)
            self.aif = self.calcAIF(60); % Gives initial shift of 60 seconds
            self.aifType = val;
        end
    end
end
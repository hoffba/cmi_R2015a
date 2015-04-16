% RegClass function
function stat = selectTform(self,varargin)
% Select current schedule step (from GUI listbox)
% Also updates GUI objects with info from ElxClass
% Syntax:
%       selectTform(#)

stat = false;

% Check for valid inputs:
if (nargin==2) && ismember(varargin{1},0:length(self.elxObj.Schedule))
    self.ind = varargin{1};
    set(self.h.listbox_Tforms,'Value',self.ind);
    stat = true;
elseif (nargin==3) && ishghandle(varargin{1}) ...
            && strcmp(get(varargin{1},'Tag'),'listbox_Tforms')
    self.ind = get(varargin{1},'Value');
    stat = true;
end

if stat 
    h = [self.h.button_SAVE,...
         self.h.button_removeTform,...
         self.h.button_AllParameters,...
         self.h.popup_Sampler,...
         self.h.popup_Metric,...
         self.h.popup_interp,...
         self.h.checkbox_saveIntermediates,...
         self.h.table_schedule,...
         self.h.edit_nres,...
         self.h.edit_defVal,...
         self.h.button_Start];
    wh = [self.h.checkbox_BEpenalty,...
          self.h.edit_finalGridX,...
          self.h.edit_finalGridY,...
          self.h.edit_finalGridZ,...
          self.h.checkbox_def,...
          self.h.checkbox_jac,...
          self.h.checkbox_jacmat];
    if isempty(self.elxObj.Schedule) || (self.ind==0)
        % Update GUI objects to disable unavailable options
        set([h,wh],'Enable','off');
    elseif (self.ind>0)
        if strcmp(get(self.h.button_SAVE,'Enable'),'off')
            % First step initaited :
            %   Need to enable relevant GUI objects
             set(h,'Enable','on');
        end
    % Set other GUI properties based on new parameters:
        p = self.elxObj.Schedule{self.ind};
    % Gray out Scales button if Transform = Warp
        str = 'off';
        if isfield(p,'Scales')
            str = 'on';
        end
        set(self.h.button_scales,'Enable',str);
    % # Resolutions:
        val = '';
        if self.ind>0
            val = num2str(self.elxObj.getPar(self.ind,'NumberOfResolutions'));
        end
        set(self.h.edit_nres,'String',val);
    % Default Value:
        val = 0;
        if isfield(p,'DefaultPixelValue')
            val = p.DefaultPixelValue;
        end
        self.h.edit_defVal.String = num2str(val);
    % ImageSampler:
        i = find(strcmp(p.ImageSampler,get(self.h.popup_Sampler,'String')));
        set(self.h.popup_Sampler,'Value',i);
    % Metric(s):
        if iscell(p.Metric)
            % Assume:   Metric{1} = main cost function
            %           Metric{2} = bending energy
            str = p.Metric{1};
            set(self.h.checkbox_BEpenalty,'Value',1);
        else
            set(self.h.checkbox_BEpenalty,'Value',0);
            str = p.Metric;
        end
        i = find(strcmp(str,get(self.h.popup_Metric,'UserData')));
        set(self.h.popup_Metric,'Value',i);
    % ResampleInterpolator:
        i = find(strcmp(p.ResampleInterpolator,get(self.h.popup_interp,'UserData')));
        set(self.h.popup_interp,'Value',i);
    % Warping grid:
        if strcmp(p.Transform,'BSplineTransform')
            set(wh,'Enable','on');
            self.h.edit_finalGridX.String = num2str(p.FinalGridSpacingInVoxels(1));
            self.h.edit_finalGridY.String = num2str(p.FinalGridSpacingInVoxels(2));
            self.h.edit_finalGridZ.String = num2str(p.FinalGridSpacingInVoxels(3));
        else
            set(wh,'Enable','off');
        end
    % Bending Energy Penalty Weight
        if isfield(p,'Metric1Weight')
            estr = 'on';
            str = num2str(p.Metric1Weight);
        else
            estr = 'off';
            str = '';
        end
        set(self.h.edit_BEpenaltyAmt,'Enable',estr,'String',str);
    % Save intermediate images:
        if isfield(p,'WriteResultImageAfterEachResolution') ...
                && strcmp(p.WriteResultImageAfterEachResolution,'true')
            i = 1;
        else
            i = 0;
        end
        set(self.h.checkbox_saveIntermediates,'Value',i);
    % UITable of options:
        n = p.NumberOfResolutions;
        % Concatenate values into table matrix:
        cstr = [reshape(p.FixedImagePyramidSchedule,3,[]);...
                reshape(p.MovingImagePyramidSchedule,3,[]);...
                p.MaximumNumberOfIterations;...
                p.SP_A;...
                p.SP_a;...
                p.NumberOfSpatialSamples];
        i = nan(3,n);
        if strcmp(p.Transform,'BSplineTransform') && isfield(p,'GridSpacingSchedule')
            i = reshape(p.GridSpacingSchedule,3,[]);
        end
        cstr = cellfun(@num2str,num2cell([cstr;i]),'UniformOutput',false);
        set(self.h.table_schedule,'Data',cstr);
    end
end

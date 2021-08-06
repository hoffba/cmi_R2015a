% RegClass function
function UDschedule(self,hObject,edata)
% Handle Callbacks from Elastix Schedule-related GUI objects

if (nargin==1)
    % Re-set Schedule to empty
    self.elxObj.rmStep;
    if self.guicheck
        set(self.h.listbox_Tforms,'String',{},'Enable','off','Value',1);
    end
elseif ishghandle(hObject)
    C = {};
    tag = get(hObject,'Tag');
    switch tag
        case 'popup_Sampler'
            str = get(hObject,'String');
            C = {'ImageSampler',str{get(hObject,'Value')}};
        case 'popup_Metric'
            str = get(hObject,'UserData');
            C = {'Metric',str{get(hObject,'Value')}};
        case 'popup_interp'
            str = get(hObject,'UserData');
            C = {'ResampleInterpolator',str{get(hObject,'Value')}};
        case 'edit_BEpenaltyAmt'
            C = {'TransformBendingEnergy',str2double(get(hObject,'String'))};
        case 'checkbox_BEpenalty'
            C = {'TransformBendingEnergy',get(hObject,'Value')};
        case 'button_removeTform'
            str = get(self.h.listbox_Tforms,'String');
            val = get(self.h.listbox_Tforms,'Value');
            stat = self.elxObj.rmStep(val);
            if stat
                str(val) = [];
                set(self.h.listbox_Tforms,'String',str);
                self.selectTform(length(self.elxObj.Schedule));
            end
        case 'button_addTform'
            self.addElxStep;
        case {'edit_finalGridX','edit_finalGridY','edit_finalGridZ'}
            i = find(strcmp(tag(end),{'X','Y','Z'}),1);
            val = str2double(get(hObject,'String'));
            n = self.elxObj.getPar(self.ind,'FinalGridSpacingInVoxels');
            if ~isnan(val)
                n(i) = val;
                C = {'FinalGridSpacingInVoxels',n};
            end
        case 'edit_nres'
            C = {'NumberOfResolutions',round(str2double(get(hObject,'String')))};
        case 'table_schedule'
            val = str2double(edata.NewData);
            if ~isnan(val)
                fldn = get(hObject,'UserData');
                fldn = fldn{edata.Indices(1)};
                vec = self.elxObj.Schedule{self.ind}.(fldn);
                switch fldn
                    case 'FixedImagePyramidSchedule'
                        i = 3*(edata.Indices(2)-1)+edata.Indices(1);
                    case 'MovingImagePyramidSchedule'
                        i = 3*(edata.Indices(2)-1)+edata.Indices(1)-3;
                    case 'GridSpacingSchedule'
                        i = 3*(edata.Indices(2)-1)+edata.Indices(1)-10;
                    otherwise
                        i = edata.Indices(2);
                end
                vec(i) = val;
                C = {fldn,vec};
            end
        case 'checkbox_saveIntermediates'
            str = 'false';
            if get(hObject,'Value')
                str = 'true';
            end
            C = {'WriteTransformParametersEachResolution',str,...
                 'WriteResultImageAfterEachResolution',str};
        case 'button_scales'
            val = self.elxObj.getPar(self.ind,'Scales');
            if ~isempty(val)
                switch length(val)
                    case 3
                        str = {'t_x','t_y','t_z'};
                    case 6
                        str = {'R_x','R_y','R_z',...
                               't_x','t_y','t_z'};
                    case 7
                        str = {'q_1','q_2','q_3',...
                               't_x','t_y','t_z','S'};
                    case 12
                        str = {'a_11','a_12','a_13',...
                               'a_21','a_22','a_23',...
                               'a_31','a_32','a_33',...
                               't_x','t_y','t_z'};
                end
                i = cellfun(@str2double,inputdlg(str,'Scales',[1,10],...
                            cellfun(@num2str,num2cell(val),...
                                    'UniformOutput',false)))';
                if ~isempty(i) && ~any(isinf(i))
                    if any(isnan(i))
                        str = 'true';
                    else
                        str = 'false';
                    end
                    C = {'Scales',i,'AutomaticScalesEstimation',str};
                end
            end
        case 'edit_defVal'
            C = {'DefaultPixelValue',str2double(get(hObject,'String'))};
        case 'checkbox_jac'
            self.jac = logical(get(hObject,'Value'));
        case 'checkbox_jacmat'
            self.jacmat = logical(get(hObject,'Value'));
        case 'checkbox_def'
            self.def = logical(get(hObject,'Value'));
        case 'checkbox_Tvoi'
            self.Tvoi = logical(get(hObject,'Value'));
        case 'checkbox_Tsurf'
            self.Tsurf = logical(get(hObject,'Value'));
        case 'checkbox_wait'
            self.waitchk = logical(hObject.Value);
        otherwise
            warning(['Unknown tag:',tag])
    end
    
    % Set ElxClass properties - Name/Value pairs in C{}
    if ~isempty(C)
        self.setElxPar(self.ind,C{:});
    end
elseif (nargin==3) && ischar(hObject)
    switch hObject
        case 'jac'
            self.jac = logical(edata);
        case 'jacmat'
            self.jacmat = logical(edata);
        case 'def'
            self.def = logical(edata);
        case 'waitchk'
            self.waitchk = logical(edata);
        case 'Tvoi'
            self.Tvoi = logical(edata);
        case 'Tsurf'
            self.Tsurf = logical(edata);
        otherwise
            warning(['Unknown tag:',tag])
    end
else
    warning('RegClass.UDschedule() : Invalid inputs');
end

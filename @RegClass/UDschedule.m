% RegClass function
function UDschedule(self,hObject,edata)
% Handle Callbacks from Elastix Schedule-related GUI objects

if ishandle(hObject)
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
            
            i = str2double(get(hObject,'String'));
            if isnan(i)
                set(hObject,'String',...
                    num2str(self.elxObj.getPar(self.ind,'Metric1Weight')));
            else
                C = {'Metric1Weight',i};
            end
            
        case 'checkbox_BEpenalty'
            
            str = self.elxObj.Schedule{self.ind}.Metric;
            if get(hObject,'Value')
                C = {'Metric',{str,'TransformBendingEnergyPenalty'},...
                     'Metric0Weight',1,'Metric1Weight',1,...
                     'UseRelativeWeights','false',...
                     'Registration','MultiMetricMultiResolutionRegistration'};
                str = {'1','on'};
            else
                C = {'Metric',str{1},...
                     'Metric0Weight',[],'Metric1Weight',[],...
                     'UseRelativeWeights',[],...
                     'Registration','MultiResolutionRegistration'};
                str = {'','off'};
            end
            set(self.h.edit_BEpenaltyAmt,'String',str{1},'Enable',str{2});
            
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
            
            opts = {'Translation','Euler','Similarity','Affine','Warp'};
            [sel,ok] = listdlg('ListString',opts,...
                               'SelectionMode','single',...
                               'Name','Transform Type');
            if ok
                answer = opts{sel};
                stat = self.elxObj.addStep(answer);
                if stat
                    if strcmp(get(self.h.listbox_Tforms,'Enable'),'off')
                        set(self.h.listbox_Tforms,'Enable','on');
                    end
                    str = [get(self.h.listbox_Tforms,'String');{answer}];
                    set(self.h.listbox_Tforms,'String',str);
                    self.selectTform(length(str));
                end
            end
            
        case {'edit_finalGridX','edit_finalGridY','edit_finalGridZ'}
            
            i = find(strcmp(tag(end),{'X','Y','Z'}));
            val = str2double(get(hObject,'String'));
            n = self.elxObj.getPar(self.ind,'FinalGridSpacingInVoxels');
            if isnan(val)
                set(hObject,'String',num2str(n(i)));
            else
                n(i) = val;
                C = {'FinalGridSpacingInVoxels',n};
            end
            
        case 'edit_nres'
            
            i = round(str2double(get(hObject,'String')));
            n = self.elxObj.getPar(self.ind,'NumberOfResolutions');
            if ~isnan(i) && (i~=n) && (i>0) && (i<6) % max of 5 resolutions for now
                % Get parameter values from table and add/remove columns from front:
                rdata = get(self.h.table_schedule,'Data');
                if i<n
                    % Remove from front
                    rdata(:,1:(n-i)) = [];
                else
                    % Add to front
                    rdata = [repmat(rdata(:,1),1,i-n),rdata];
                end
                set(self.h.table_schedule,'Data',rdata);
                rdata = cellfun(@str2double,rdata);
                % Update table-relevant Schedule options:
                self.elxObj.setPar(self.ind,...
                   'NumberOfResolutions',i,...
                   'FixedImagePyramidSchedule',reshape(rdata(1:3,:),[1,i*3]),...
                   'MovingImagePyramidSchedule',reshape(rdata(4:6,:),[1,i*3]),...
                   'MaximumNumberOfIterations',rdata(7,:),...
                   'SP_A',rdata(8,:),...
                   'SP_a',rdata(9,:),...
                   'NumberOfSpatialSamples',rdata(10,:));
               if strcmp(self.elxObj.getPar(self.ind,'Transform'),'BSplineTransform')
                   self.elxObj.setPar(self.ind,'GridSpacingSchedule',...
                       reshape(rdata(11:13,:),[1,i*3]));
               end
            end
            set(hObject,'String',...
                num2str(self.elxObj.getPar(self.ind,'NumberOfResolutions')));
            
        case 'table_schedule'
            
            val = str2double(edata.NewData);
            if isnan(val) % Set value back to original
                cstr = get(hObject,'Data');
                cstr{edata.Indices(1),edata.Indices(1)} = edata.PreviousData;
                set(hObject,'Data',cstr);
            else
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
            
            val = self.elxObj.Schedule{self.ind}.Scales;
            if ~isempty(val)
                switch self.elxObj.getPar(self.ind,'Transform')
                    case 'TranslationTransform'
                        str = {'t_x','t_y','t_z'};
                    case 'EulerTransform'
                        str = {'R_x','R_y','R_z',...
                               't_x','t_y','t_z'};
                    case 'SimilarityTransform'
                        str = {'q_1','q_2','q_3',...
                               't_x','t_y','t_z','S'};
                    case 'AffineTransform'
                        str = {'a_11','a_12','a_13',...
                               'a_21','a_22','a_23',...
                               'a_31','a_32','a_33',...
                               't_x','t_y','t_z'};
                end
                i = cellfun(@str2double,inputdlg(str,'Scales',[1,10],...
                            cellfun(@num2str,num2cell(val),...
                                    'UniformOutput',false)))';
                if ~isempty(i) && ~any(isnan(i) | isinf(i))
                    C = {'Scales',i};
                end
            end
            
        case 'edit_defVal'
            
            val = str2double(get(hObject,'String'));
            if isnan(val)
                set(hObject,'String',num2str(self.elxObj.getPar(self.ind,'DefaultPixelValue')));
            else
                C = {'DefaultPixelValue',val};
            end
            
        case 'checkbox_tformVOI'
            
            self.tVOI = logical(get(hObject,'Value'));
            
        case 'checkbox_jac'
            
            self.jac = logical(get(hObject,'Value'));
            
        case 'checkbox_jacmat'
            
            self.jacmat = logical(get(hObject,'Value'));
            
        case 'checkbox_def'
            
            self.def = logical(get(hObject,'Value'));
            
        otherwise
            warning(['Unknown tag:',tag])
    end
    
    % Set ElxClass properties - Name/Value pairs in C{}
    if ~isempty(C)
        self.elxObj.setPar(self.ind,C{:});
    end
else
    warning('RegClass.UDschedule() : Invalid inputs');
end

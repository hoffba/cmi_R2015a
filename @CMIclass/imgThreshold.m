% CMIcalss function
% Threshold Image
function imgThreshold(self,hObject,~)
% GUI callback to set image thresholds and update displayed image/plot

if self.img.check
    go = true;
    cvec = self.vec;
    tval = self.img.thresh(self.vec,:)';
    vlims = self.getProp('ValLim');
    switch hObject.Tag
        case 'checkbox_thMin'
            if hObject.Value
                tval(1) = self.img.valExt(cvec,1);
                set(self.h.slider_thMin,'Enable','on','Min',vlims(1),...
                    'Max',min(vlims(2),tval(2)),'Value',vlims(1));
                set(self.h.edit_thMin,'Enable','on','String',num2str(vlims(1)));
            else
                tval(1) = -inf;
                set([self.h.slider_thMin,self.h.edit_thMin],'Enable','off');
            end
        case 'checkbox_thMax'
            if hObject.Value
                tval(2) = self.img.valExt(self.vec,2);
                set(self.h.slider_thMax,'Enable','on','Min',max(vlims(1),tval(1)),...
                    'Max',vlims(2),'Value',vlims(2));
                set(self.h.edit_thMax,'Enable','on','String',num2str(vlims(2)));
            else
                tval(2) = inf;
                set([self.h.slider_thMax,self.h.edit_thMax],'Enable','off');
            end
        case 'slider_thMin'
            tval(1) = hObject.Value;
            set(self.h.slider_thMax,'Min',tval(1));
            set(self.h.edit_thMin,'String',num2str(tval(1)));
        case 'slider_thMax'
            tval(2) = hObject.Value;
            set(self.h.slider_thMin,'Max',tval(2));
            set(self.h.edit_thMax,'String',num2str(tval(2)));
        case 'edit_thMin'
            a = str2double(hObject.Value);
            if isnan(a)
                go = false;
                set(hObject,'String',num2str(tval(1)));
            else
                tval(1) = a;
                set(self.h.slider_thMax,'Min',tval(1));
            end
        case 'edit_thMax'
            a = str2double(hObject.Value);
            if isnan(a)
                go = false;
                set(hObject,'String',num2str(tval(2)));
            else
                tval(2) = a;
                set(self.h.slider_thMin,'Max',tval(2));
            end
        case 'image_thresh'
            defAns = cellstr(num2str(tval'))';
            answer = inputdlg({'Threshold MIN value:','Threshold MAX value:'},...
                'Set image threshold',1,defAns);
            tval = str2double(answer);
        otherwise
            go = false;
    end

    if go
        if self.applyallcheck
            cvec = 1:self.img.dims(4);
            tval = tval * ones(1,length(cvec));
        end
        self.img.threshold(cvec,tval);
        self.dispUDslice;
    end
end
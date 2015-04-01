% CMIcalss function
% Delete Image
function imgDelete(self,~,~)
nvec = self.img.dims(4);
if (nvec > 1) % don't delete if only 1 image
    inds = listdlg('ListString',self.img.labels,...
        'PromptString','Image vector(s) to remove:');
    % Do nothing unless all input values are valid
    if ~isempty(inds)
        self.img.imgDelete(inds);
        self.clim(inds,:) = [];
        if any(self.vec == inds)
            % Reset current image to first
            self.vec = 1;
        else
            % Shift vec with matrix location
            self.vec = self.vec - numel(find(inds < self.vec));
        end
        nv = self.img.dims(4);
        cvec = min(self.vec,nv);
        if nv > 1
            set(self.h.slider_4D,'Value',cvec,'Max',nv,...
                                 'SliderStep',[1,1]/(nv-1),'Enable','on');
        else
            set(self.h.slider_4D,'Enable','off');
        end
        set(self.h.text_nv,'String',num2str(nv));
        set(self.h.edit_vec,'String',num2str(cvec));
        self.setVec(cvec);
    end
end
% CMIclass function
function ind = imgAppend(self,img,label)

if (nargin>1) && ~isempty(img)
    [d(1),d(2),d(3),d(4)] = size(img);
    if all(d(1:3)==self.img.dims(1:3))
        if (nargin<3) || ~iscellstr(label)
            label = strcat({'Append-'},num2str((1:d(4))','%02u'))';
        end
        ind = self.img.imgAppend(img,label);
        nvec = size(self.clim,1)+1 : self.img.dims(4);
        self.clim(nvec,:) = [squeeze(min(min(min(self.img.mat(:,:,:,nvec))))),...
                             squeeze(max(max(max(self.img.mat(:,:,:,nvec)))))];
        d = self.img.dims(4);
        estr = 'off';
        if d>1
            estr = 'on';
        end
        set(self.h.text_nv,'String',num2str(d));
        set(self.h.slider_4D,'Max',d,'SliderStep',[1,1]/(d-1),'Enable',estr);
        set(self.h.edit_vec,'Enable',estr);
        self.setVec(d);
    end
end


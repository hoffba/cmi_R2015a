% ImageClass function
% Save image
function fname = saveImg(self,tvec,fname,cvec)
if self.check
    status = true;
    if isempty(fname)
        [~,fname,~] = fileparts(self.name);
        %fname = self.name;
    end
    % Save full 4D or just current image?
    if ((nargin<2) || isempty(tvec))
        if self.dims(4)>1
            answer = questdlg('Save full 4D or only current image?','Save Image',...
                'All','Current','Cancel','All');
            switch answer
                case 'All'
                    tvec = 1:self.dims(4);
                case 'Current'
                    tvec = cvec;
                case 'Cancel'
                    status = false;
            end
        else
            tvec = cvec;
        end
    end
    % If OK, save image to file
    if status
        fov = self.dims(1:3).*self.voxsz;
        if (nargin<3) % fname was not input
            fname = self.name;
        end
        status = cmi_save(0,self.mat(:,:,:,tvec),self.labels(tvec),fov,fname);
        if status
            fname = status;
            [self.dir,self.name] = fileparts(status);
        else
            fname = '';
        end
    else
        fname = '';
    end
else
    fname = '';
end
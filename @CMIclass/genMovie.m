% CMIclass function
% Create Movie of current image
function genMovie(self,~,~)
if self.img.check && self.guicheck
    F = self.grabFrames;
    if ~isempty(F)
        % Write each padded frame to the video file
        frate = str2double(inputdlg({'Desired Frame Rate:'},'Frame Rate',1,{'5'}));
        if ~isnan(frate)
            [fname,path] = uiputfile('*.avi','Save Movie As');
            fname = fullfile(path,fname);
            writerObj = VideoWriter(fname);
            writerObj.FrameRate = frate;
            open(writerObj);
            [nx,ny,~] = size(F(1).cdata);
            timg = zeros(nx,ny,3);
            for i = 1:length(F)
                timg = 0*timg;
                tF = F(i);
                tnx = size(tF.cdata,1);
                tny = size(tF.cdata,2);
                padx = round((nx-tnx)/2);
                pady = round((ny-tny)/2);
                timg((padx+1):(padx+tnx),(pady+1):(pady+tny),:) = tF.cdata;
                tF.cdata = uint8(timg);
                writeVideo(writerObj,tF);
            end
        close(writerObj);
        end
    end
end
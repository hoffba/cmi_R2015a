% RegClass function
function stat = oneshot(self,~,~)
% stat = oneshot;
%           Uses pre-set GUI and RegClass properties (e.g. RegClass.points)

gochk = false;
if ~(self.cmiObj(1).img.check && self.cmiObj(2).img.check)
    error('Need to load images before transforming.')
elseif get(self.h.checkbox_useExistingT,'Value')
    gochk = true;
else
    % Calculate new Tform based on existing points:
    inC = {'DefaultPixelValue',self.T0defVal};
    if get(self.h.popup_Transforms,'Value') == 6 % Warp
%         np = min(size(self.points{1},1),size(self.points{2},1));
%         off1 = self.cmiObj(1).img.voxsz.*self.cmiObj(1).img.dims(1:3)/2;
%         off2 = self.cmiObj(2).img.voxsz.*self.cmiObj(2).img.dims(1:3)/2;
%         p1 = ((self.points{1}(1:np,:)-repmat(off1([2,1,3]),np,1)))';
%         p2 = ((self.points{2}(1:np,:)-repmat(off2([2,1,3]),np,1)))';
%         M = reshape(p2',1,[]);
%         fM = reshape(p1',1,[]);
        M = reshape(self.ctrpts(2)',1,[]);
        fM = reshape(self.ctrpts(1)',1,[]);
        inC = [inC,{'FixedImageLandmarks',fM}];
    else
        M = self.pts2M; % Use existing points, so empty input
        M = [reshape(M(1:3,1:3)',1,[]),M(1:3,4)'];
    end
    if ~isempty(M)
        % Set elxObj Initial Transform properties:
        self.elxObj.setTx0(M,...
                self.cmiObj(1).img.voxsz([2,1,3]),...
                self.cmiObj(1).img.dims([2,1,3]),...
                inC{:});
        set(self.h.checkbox_useExistingT,'Enable','on','Value',1);
        gochk = true;
    end
end

if gochk

    % Check for ouput directory:
    if isempty(self.odir)
        self.setOdir;
    end
    stat = exist(self.odir,'dir');

    % Call transformix:
    if stat
        hw = waitbar(0,'Saving temporary files for Transformix ...');
        tpfname = self.elxObj.saveTx0(fullfile(self.odir,...
            ['TransformParameters-snap-',datestr(now,'yyyymmdd'),'.txt']));
        fname = fullfile(self.odir,'elxtemp-in.nii');
        stat = saveNIFTI(fname,self.cmiObj(2).img.mat(:,:,:,self.cmiObj(2).vec),...
                       [],self.cmiObj(2).getProp('fov'),self.cmiObj(2).img.orient);
        if stat
            waitbar(0.5,hw,'Calling Transformix ...');
            
            % Set Reference image properties:
            self.elxObj.setTx0par('Size',self.cmiObj(1).img.dims([2,1,3]),...
                                  'Spacing',self.cmiObj(1).img.voxsz([2,1,3]));

            % Generate system command:
            cmdstr = self.elxObj.sysCmd(self.odir,'title','Transformix',...
                'tp',tpfname,'in',fname,'wait',true);
            % Save command to file for trouble-shooting:
            fid = fopen(fullfile(self.odir,'transformixCMD.txt'),'w');
            if fid>2
                fprintf(fid,'%s',cmdstr);
                fclose(fid);
            end
            stat = ~system(cmdstr);%csh''&']);

            % Append transformed image to Reference CMIobj
            if stat
                waitbar(0.9,hw,'Appending transformed image to Reference ...');
                self.cmiObj(1).loadImg(true,fullfile(self.odir,'result.nii'));
            end
        end
        delete(hw);
    end
end
    
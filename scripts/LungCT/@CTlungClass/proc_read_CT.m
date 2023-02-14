function [stat,swap_flag] = proc_read_CT(self,dat_str)

swap_flag = false;
if nargin==1 || isempty(dat_str)
    dat_str = {'ct_exp','ct_ins'};
elseif ~all(ismember(dat_str,{'ct_exp','ct_ins'}))
    return;
end

% Determine file names:
fn_save = self.getFileName(dat_str);

N = numel(dat_str);
img_chk = false(1,N);
for i = 1:N
    if exist(fn_save{i},'file')
        writelog(self.fn_log,'Loading %s from file.\n',self.data.(dat_str{i}).tag);
        self.dat.(dat_str{i}).info = niftiinfo(fn_save{i});
        self.dat.(dat_str{i}).mat = niftiread(fn_save{i});
        img_chk(i) = true;
    else
        writelog(self.fn_log,'Loading %s from DICOM.\n  Path: ',self.dat.(dat_str{i}).tag);
        if isempty(self.dat.(dat_str{i}).dcmpath)
            warning('Valid DICOM path not found\n');
        else
            writelog(self.fn_log,'%s\n',self.dat.(dat_str{i}).dcmpath);
            % Load DICOM image
            [self.dat.(dat_str{i}).mat,~,fov,orient,~] = readDICOM(dcmpath);
            d = size(self.dat.(dat_str{i}).mat);
            voxsz = fov./d;

            % Correct image orientation by bone threshold
            if self.orient_flag
                prop = regionprops(max(self.dat.(dat_str{i}).mat(:,:,round(d(3)/2):end)>800,[],3),'Orientation','Area');
                if mod(round(prop([prop.Area]==max([prop.Area])).Orientation/90),2)
                    writeLog(self.fn_log,'   Permuting image\n');
                    ind = [2,1,3];
                    self.dat.(dat_str{i}).mat = permute(self.dat.(dat_str{i}).mat,ind);
                    voxsz = voxsz(ind);
                    orient = orient([ind,4],[ind,4]);
                end
            end

            % Setup NIFTI info
            self.dat.(dat_str{i}).info = init_niftiinfo(self.dat.(dat_str{i}).tag,voxsz,'int16',d);
            self.dat.(dat_str{i}).info.Transform = affine3d((diag([-1 -1 1 1])*orient)');
            
            img_chk(i) = true;
        end
    end
end

% Check EXP/INSP:
if self.swap_flag && N==2 && all(img_chk)
    
    % Quick segmentation for Exp/Ins Identification:
    ref_seg = getRespiratoryOrgans(medfilt2_3D(self.dat.ct_ref.mat));
    hom_seg = getRespiratoryOrgans(medfilt2_3D(self.dat.ct_hom.mat));

    voxvol = [prod(self.dat.ct_ref.info.PixelDimensions) , prod(self.dat.ct_hom.info.PixelDimensions)];
    
    if (nnz(ref_seg)*voxvol(1)) > (nnz(hom_seg)*voxvol(2))
        swap_flag = true;
        writeLog(self.fn_log,'Swapping INS/EXP due to lung volume\n');
        timg = self.dat.ct_ref;
        self.dat.ct_ref = self.dat.ct_hom;
        self.dat.ct_hom = timg;
        
        % Don't swap the tags
        self.dat.ct_hom.tag = self.dat.ct_ref.tag;
        self.dat.ct_ref.tag = timg.tag;
    end
end

% Save images as nii
for i = 1:N
    writeLog(self.fn_log,'Saving image: %s\n',fn_save{i});
    niftiwrite(self.dat.(dat_str{i}).mat,fn_save{i}(1:end-3),self.dat.(dat_str{i}).info,'Compressed',true);
end

stat = all(img_chk);

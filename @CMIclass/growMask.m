% CMIclass function
% segmentation algorithms for mask region growing
function growMask(self,~,~)
if self.img.check
    opts = {'Deep Learning','Points-Threshold','VOI-Threshold','Ends-in','Lung Segment'};
    [sel,ok] = listdlg('PromptString','Select a mask growing option:','ListString',opts);
    if ok
        % Threshold for region growing is set by user before running this algorithm
        switch sel
            case 1 % Deep Learning
                tmask = DL_lung_segmetation(self.img.mat(:,:,:,self.vec));
                self.img.mask.merge('replace',logical(tmask));
                self.imgAppend(tmask,{'Segmentation'});
            case 2 % Points-Threshold
                stat = true;
                if self.img.mask.check
                    stat = strcmp(questdlg('This will delete your current mask. Continue?'),'Yes');
                end
                if stat
                    [x,y] = getpts(self.haxes);
                    tdims = self.img.dims(1:3);
                    tdims(self.orient) = [];
                    tvsz = self.img.voxsz;
                    tvsz(self.orient) = [];
                    x = round(x/tvsz(2));
                    y = round(y/tvsz(1));
                    tmask = false(tdims);
                    for i = 1:length(x)
                        tmask(y(i),x(i)) = true;
                    end
                    self.img.mask.clear;
                    self.img.mask.setSlice(tmask,self.orient,1,self.slc(self.orient));
                    self.img.growVOI(self.vec);
                end
            case 3 % VOI-Threshold
                if any(self.img.mask.mat(:))
                    self.img.growVOI(self.vec)
                else
                    error('Need to draw a mask first!')
                end
            case 4 % Ends-in
                self.img.growVOI(self.vec,1);
            case 5 % Automated Lung Segmentation
                lmask = self.img.segmentLung(self.vec);
                if ~islogical(lmask)
                    self.imgAppend(lmask,{'Segmentation'});
                end
        end
        self.dispUDmask;
    end
end
% CMIclass function
% segmentation algorithms for mask region growing
function growMask(self,~,~)
if self.img.check
    opts = {'Points-Threshold','VOI-Threshold','Ends-in','Lung Segment'};
    [sel,ok] = listdlg('PromptString','Select a mask growing option:','ListString',opts);
    if ok
        % Threshold for region growing is set by user before running this algorithm
        switch sel
            case 1 % Points-Threshold
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
            case 2 % VOI-Threshold
                if any(self.img.mask.mat(:))
                    self.img.growVOI(self.vec)
                else
                    error('Need to draw a mask first!')
                end
            case 3 % Ends-in
                self.img.growVOI(self.vec,1);
            case 4 % Automated Lung Segmentation
                self.img.segmentLung(self.vec);
        end
        self.dispUDmask;
    end
end
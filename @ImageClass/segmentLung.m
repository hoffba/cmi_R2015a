% ImageClass function
function segmentLung(self,vec,humanchk)
lmask = [];
if (nargin<3)
    humanchk = questdlg('Choose Data Type:','Lung Segmentation',...
                            'Mouse','Human','Mouse');
    humanchk = strcmp(humanchk,'Human');
end

if humanchk % human
        lmask = segLungHuman(self.mat(:,:,:,vec),self.mask.mat);
else % mouse
    if (nargin==2)
        answer = inputdlg({'Trachea threshold:','Lung threshold:',...
                        'Noise Filter:','Dilation:','Erosion:'},...
                        'Lung Segmentation',1,{'150','700','5','3','2'});
        answer = str2double(answer);
    else
        answer = [150,700,5,3,2];
    end
    if ~isempty(answer) && ~any(isnan(answer))
        lmask = segLungMouse(self.mat(:,:,:,vec),...
                             struct('Tt',{answer(1)},...
                                    'Tl',{answer(2)},...
                                    'r' ,{answer(3)},...
                                    'me',{answer(4)},...
                                    'md',{answer(5)}));
    end
end
if ~isempty(lmask)
    self.mask.merge('Replace',lmask);
end

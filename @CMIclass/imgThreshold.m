% CMIcalss function
% Threshold Image
function imgThreshold(self,~,~)
if self.img.check
    defAns = cellstr(num2str(self.img.thresh(self.vec,:)'))';
    answer = inputdlg({'Threshold MIN value:','Threshold MAX value:'},...
        'Set image threshold',1,defAns);
    tval = str2double(answer);
    if (length(tval) == 2)
        if tval(2) > tval(1)
            tvec = self.vec;
            if self.applyallcheck
                tvec = 1:self.img.dims(4);
                tval = tval * ones(1,length(tvec));
            end
            self.img.threshold(tvec,tval);
            self.dispUDslice;
        end
    end
end
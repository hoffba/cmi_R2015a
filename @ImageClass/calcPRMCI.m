% ImageClass function
function out = calcPRMCI(self,vec)
if self.check && self.mask.check
    tvec = self.prm.dvec;
    tvec(tvec==0) = vec;
    if length(tvec)>1
        
        % Grab masked image values:
        Pixels = self.getMaskVals(tvec(1:2));
        
        % Use PRM cutoffs:
        tcut = self.prm.cutoff; tcut(tcut(:,1)==0,1) = vec;
        if ~isempty(tcut)
            ind = (Pixels(:,tvec(1))<tcut(tvec(1)==tcut(:,1),2)) ...
                | (Pixels(:,tvec(1))>tcut(tvec(1)==tcut(:,1),3)) ...
                | (Pixels(:,tvec(2))<tcut(tvec(2)==tcut(:,1),2)) ...
                | (Pixels(:,tvec(2))>tcut(tvec(2)==tcut(:,1),3));
            Pixels(ind,:) = [];
        end
        [slope,YInt,CI95] = PRM_threshold2(Pixels);
        out = struct('slope',slope,'Yint',YInt,'CI95',CI95);
        assignin('base','prmCI',out);
    end
end
% ImageClass function
function stat = selectPRMscatter(self,nres)

if nargin==1
    nres = 200;
end
stat = self.prm.scatterSelect( self.mask.mat,...
                               self.getMaskVals(self.prm.dvec(1:2)),...
                               nres );

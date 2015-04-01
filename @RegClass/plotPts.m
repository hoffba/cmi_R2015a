% RegClass function
function plotPts(self,varargin)
% Plot user-selected points on CMIcalss figure

h = [];
if nargin==2
    % Manual callback
    % varargin = { image#(1/2,Ref/Hom) }
    i = varargin{1};
    h = self.cmiObj(i).haxes;
elseif (nargin==3)
    % GUI listener callback
    h = varargin{2}.AffectedObject.haxes;
    i = find(cellfun(@(x)~(isempty(x)||isempty(h))&&(x==h),{self.cmiObj(:).haxes}),1);
end

if ~isempty(h)
    % Convert points back to corner-centric geometry
    p = self.points{i};
    if ~isempty(p)
        ind = round(0.5 + p(:,self.cmiObj(i).orient)/...
                            self.cmiObj(i).img.voxsz(self.cmiObj(i).orient))...
              == self.cmiObj(i).slc(self.cmiObj(i).orient);
        p = p(ind,(1:3)~=self.cmiObj(i).orient);
    end
    if isempty(p)
        if ~isnan(self.hpts(i)) && ishghandle(self.hpts(i))
            delete(self.hpts(i));
            self.hpts(i) = nan;
        end
    else
        if isnan(self.hpts(i)) || ~ishghandle(self.hpts(i))
            hold(h,'on');
            self.hpts(i) = plot(h,p(:,2),p(:,1),'*g');
            hold(h,'off');
        else
            set(self.hpts(i),'XData',p(:,2),'YData',p(:,1));
        end
    end
end
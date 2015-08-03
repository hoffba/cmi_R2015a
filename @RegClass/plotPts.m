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
    % POINTS ARE IN XYZ
    p = self.points{i};
    if ~isempty(p)
        
        ornt = [2,1,3];
        ornt = ornt(self.cmiObj(i).orient);
        fov = self.cmiObj(i).getProp('fov');
        voxsz = self.cmiObj(i).img.voxsz;
        
        % Convert slice coordinates to matrix indices and find those on current slice:
        ind = round((p(:,ornt) + fov(ornt)/2)/voxsz(ornt) + 0.5) == self.cmiObj(i).getProp('slc');
        
        % Find points on current slice:
        ind2 = (1:3)~=ornt;
        p = bsxfun(@plus,p(ind,ind2),fov(ind2)/2);
        
    end
    if isempty(p)
        if ~isnan(self.hpts(i)) && ishghandle(self.hpts(i))
            delete(self.hpts(i));
            self.hpts(i) = nan;
        end
    else
        if isnan(self.hpts(i)) || ~ishghandle(self.hpts(i))
            hold(h,'on');
            self.hpts(i) = plot(h,p(:,1),p(:,2),'*g');
            hold(h,'off');
        else
            set(self.hpts(i),'XData',p(:,1),'YData',p(:,2));
        end
    end
end
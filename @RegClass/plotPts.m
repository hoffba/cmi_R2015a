% RegClass function
function plotPts(self,varargin)
% Plot user-selected points on CMIclass figure

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
    % Points are stored in rectilinear xyz
    p = self.points{i};
    if ~isempty(p)
        
        % Retrieve image properties:
        ornt = self.cmiObj(i).orient;
        slc = self.cmiObj(i).slc(ornt);
        voxsz = self.cmiObj(i).img.voxsz;
        
        % Find points on this slice
        xyi = [2 1 3];
        p(round(p(:,xyi(ornt))/voxsz(xyi(ornt)))~=slc,:) = [];
        xyi(ornt) = []; xyi = flip(xyi);
        p = p(:,xyi);
        
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
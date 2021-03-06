% RegClass function
function addPoint(self,hObject,edata)
% Add selected point to list
% addPoint(image#,pts)
%       image# = 1 or 2 (Ref or Hom)
%       pts = [N x 3] List of 3D coordinates to add (X,Y,Z)
% ** Remember that Elasix geometry puts (0,0,0) at center of 3D FOV!

ci = [];
if nargin==3
    if isa(edata,'matlab.graphics.eventdata.Hit')
        % Check for Ctrl+click:
        if strcmp(get(get(get(hObject,'Parent'),'Parent'),'SelectionType'),'alt')
            % GUI call
            % Determine which CMI Object it came from:
            ci = find(cellfun(@(x)~isempty(x)&&(x==hObject),...
                                {self.cmiObj(:).hiover}),1);
            tObj = self.cmiObj(ci);
            
            i = 1:3;
            i(tObj.orient) = [];
            pos = get(get(hObject,'Parent'),'CurrentPoint')+0.5;
            ind([i,tObj.orient]) = [ pos(1,2) , pos(1,1) , tObj.slc(tObj.orient) ];
            p = tObj.img.getImageCoords(ind);

        end
    elseif ismember(hObject,[1,2]) && (size(edata,2)==3)
        % Manual input points [x,y,z]
        ci = hObject;
        p = edata;
    end
    if ~(isempty(ci) || isempty(p))
        self.points{ci}(end+(1:size(p,1)),:) = p;
        self.plotPts(ci);
    end
end


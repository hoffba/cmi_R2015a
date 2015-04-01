% RegClass function
function addPoint(self,hObject,edata)
% Add selected point to list
% addPoint(image#,pts)
%       image# = 1 or 2 (Ref or Hom)
%       pts = [N x 3] List of 3D coordinates to add
% ** Remember that Elasix geometry puts (0,0,0) at center of 3D FOV!

i = [];
if nargin==3
    if isa(edata,'matlab.graphics.eventdata.Hit')
        % Check for Ctrl+click:
        if strcmp(get(get(get(hObject,'Parent'),'Parent'),'SelectionType'),'alt')
            % GUI call
            % Determine which CMI Object it came from:
            i = find(cellfun(@(x)~isempty(x)&&(x==hObject),...
                                {self.cmiObj(:).hiover}),1);
            tObj = self.cmiObj(i);
            p = get(get(hObject,'Parent'),'CurrentPoint');
            % [x,y] --> [row,column,slice] dimensions
            p = [p(1,[2,1]),(tObj.slc(tObj.orient)-0.5)*tObj.img.voxsz(tObj.orient)];
            if tObj.orient==2
                p = p([1,3,2]);
            elseif tObj.orient==1
                p = p([3,1,2]);
            end
            % Shift to image-centric geometry
%             p = p - tObj.img.dims(1:3).*tObj.img.voxsz/2;
        end
    elseif ismember(hObject,[1,2]) && (size(edata,2)==3)
        % Manual input points [x,y,z]
        i = hObject;
        p = edata;
    end
    if ~(isempty(i) || isempty(p))
        self.points{i}(end+(1:size(p,1)),:) = p;
        self.plotPts(i);
    end
end


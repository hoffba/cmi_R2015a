% RegClass function
function addPoint(self,hObject,edata)
% Add selected point to list
% addPoint(image#,pts)
%       image# = 1 or 2 (Ref or Hom)
%       pts = [N x 3] List of 3D coordinates to add (X,Y,Z)
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
            ornt = tObj.orient;
            
            % Grap selected 2D point from axes (x,y):
            p = get(get(hObject,'Parent'),'CurrentPoint');
            p = [ p(1,1:2) , (tObj.slc(ornt)-0.5)*tObj.img.voxsz(ornt) ];
            
            % Permute to original matrix coordinates
            switch ornt
                case 1
                    torder = [2,3,1];
                case 2
                    torder = [3,2,1];
                case 3
                    torder = [1,2,3];
            end
            
            % convert to image-centric coordinates
            fov = tObj.getProp('fov');
            p = p(torder) - fov([2,1,3])/2;
            
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


% RegClass function
function addPoint(self,hObject,edata)
% Add selected point to list
% addPoint(image#,pts)
%       image# = 1 or 2 (Ref or Hom)
%       pts = [N x 3] List of 3D coordinates to add (X,Y,Z)

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
            voxsz = tObj.img.voxsz(ornt);
            
            % Grab selected 2D point from axes (x,y)
            %       (use ImageClass.orient to convert to referenced spatial coordinates)
            p = get(get(hObject,'Parent'),'CurrentPoint');
            
            % Determine coordinate order
            torder = [2 1 3]; t = torder(ornt); torder(ornt) = [];
            torder = flip(torder); torder(end+1) = t;
            [~,torder] = sort(torder);

            % Permute point to xyz
            p = [ p(1,1:2) (tObj.slc(ornt)-0.5)*voxsz ];
            p = p(torder);
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


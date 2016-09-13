% ImageIO class function
% stat = write(M,info)
% stat = write(M,info,fname)

function stat = write(self,M,fname,info)

stat = false;
if (nargin>3) && ismatrix(M) && ischar(fname)
    if (nargin<5) || ~isstruct(info) || ~all(isfield(info,{'voxsz'}))
        info = struct('VoxelSize',ones(1,3),'VoxelSpacing',ones(1,3),...
                       'Orientation',[1,0,0,0,1,0,0,0,1],...
                       'AnatomicOrientation','HFS');
    end
    if (nargin<6) || ~iscellstr(labels)
        labels = cellfun(@num2str,num2cell(1:size(M,4)),'UniformOutput',false);
    end
    [d(1),d(2),d(3)] = size(M);
    stat = feval(self.getFunc(fname,1),fname,M,labels,d.*geo.voxsz,info);
end
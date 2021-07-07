% ImageClass function
function setMat(self,imat,ilabels,ivoxsz,iorient,iname)
% Manually set the matrix
d = size(imat);
if all(d>0) && (length(d)>2)
    if length(d)>4
        d = d(1:4);
        imat = reshape(imat(1:prod(d)),d);
    elseif length(d)<4
        d(4) = 1;
    end
    self.prmBaseVec = 1;
    self.mat = imat;
    self.valExt = [squeeze(min(min(min(imat,[],1),[],2),[],3)),...
                   squeeze(max(max(max(imat,[],1),[],2),[],3))];
    self.dims = d;
    self.orient = iorient;
    self.mask.setDims(d(1:3));
    self.thresh = (inf*[-1 1]'*ones(1,d(4)))';
    self.scaleM = ones(1,d(4));
    self.scaleB = zeros(1,d(4));
    self.check = true;
    if (nargin<3) || ~iscellstr(ilabels) || (length(ilabels)~=d(4))
        ilabels = strcat('vec',cellfun(@num2str,num2cell(1:d(4)),...
                                        'UniformOutput',false));
    end
    self.labels = ilabels;
    if (nargin<4) || ~isnumeric(ivoxsz) || (length(ivoxsz)~=3)
        ivoxsz = ones(1,3);
    end
    self.voxsz = ivoxsz;
    if (nargin<5) || ~ischar(iname)
        iname = 'User-Input Image';
        idir = '';
    else
        [idir,iname] = fileparts(iname);
    end
    self.name = iname;
    self.dir = idir;
    self.mask.initialize(d(1:3));
    self.prm.initialize;
end
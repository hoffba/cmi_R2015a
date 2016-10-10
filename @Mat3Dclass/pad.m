% Mat3Dclass function
function pad(self,dim,r,val)
% Pad the current image
% Inputs: 
%   dim = vector of dimensions to pad (i.e. 1 --> rows)
%   r   = pad size for each dimension
%   val = [cell array] pad value for each image
%           (numeric, 'circular', 'replicate', or 'symmetric')

if ~isempty(self.mat) && isnumeric(dim) && isnumeric(r) && iscell(val)
    d4 = self.dims(4);
    
    % May set all values the same with one input:
    nv = length(val);
    nd = length(dim);
    nr = length(r);
    if nr==1
        r = r*ones(1,nd);
    end
    if nv==1
        val = repmat(val,1,d4);
        nv = d4;
    end
    
    % Validate inputs:
    if ~all(ismember(dim,1:3))
        error('Invalid pad dimension. Must be 1:3.');
    elseif nv~=d4
        error(['Wrong number of values input for padding (',...
            num2str(nv),'), need ',num2str(d4)]);
    else
        
        % Determine new dimensions:
        padsize = zeros(1,3);
        padsize(dim) = r;
        dout = [ 2*padsize + self.dims(1:3) , d4 ];
        
        % Loop over each image:
        tmat = nan(dout);
        try
        for i = 1:d4
            tmat(:,:,:,i) = padarray(self.mat(:,:,:,i),padsize,val{i});
        end
        catch err
            disp('check')
        end
        
        % Set new properties:
        self.mat = tmat;
        self.dims = dout;
        
        % Also pad Mat3Dclass properties of this class:
        if isprop(self,'mask') && isa(self.mask,'MaskClass') && isvalid(self.mask)
            % Always pad with 0
            self.mask.pad(dim,r,{false});
        end
        if isprop(self,'prm') && isa(self.prm,'PRMclass') && isvalid(self.prm)
            % Always pad with 0
            self.prm.pad(dim,r,{0});
        end
    end
end
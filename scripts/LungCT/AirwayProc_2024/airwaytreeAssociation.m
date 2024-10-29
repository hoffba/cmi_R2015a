% Assigns quantitative values to distal nodes in the airway tree
function [V,V_label] = airwaytreeAssociation(B,N,voxsz,r,M,M_label)
% Inputs:
%   B = Branches matrix
%   N = Nodes matrix
%   voxsz = voxel dimensions
%   r = radius of window
%   M = quantitative maps to associate with airway tree (4D)
%   M_label = {Nx2} cell array of M labels (col1) and function flags (col2)
% Function flags:
%   1 = mean
%   2 = %Volume of binary matrix
%   2 = PRM% (5-color)

%% Initialization
d = size(M,1:3);
np = prod(d);
nv = size(M,4);

% Window for grabbing voxels values around a node (n x 3 matrix)
se = calcSe(voxsz,r);

% Find terminal branches
iterm = find(B(:,6)==1);
nterm = numel(iterm);

% Convert node locations to indices
Ni = round(N(:,2:4)./voxsz);

% Determine number of output values for each node
prm_val = cell(nv,1);
n = 0;
for i = 1:nv
    fcn_ind = M_label{i,2};
    if ismember(1,fcn_ind)
        n = n+1;
    end
    if ismember(2,fcn_ind)
        n = n+1;
    end
    if ismember(3,fcn_ind)
        uM = unique(M(:,:,:,i));
        uM(uM==0) = [];
        prm_val{i} = uM;
        n = n + numel(uM);
    end
end

%% Associate quantitative values to node locations using window
ii = 0;
V = [ B(iterm,1) , nan(nterm,n) ];
V_label = strings(1,n);
for i = 1:nterm

    % Determine window locations in matrix centered on distal node
    w = se + Ni(B(iterm(i),3),:);

    % Remove invalid indices
    ind = any(w<1,2) | w(:,1)>d(1) | w(:,2)>d(2) | w(:,3)>d(3);
    w(ind,:) = [];

    % Convert window subscripts to indices
    ind = sub2ind(d,w(:,1),w(:,2),w(:,3));

    % Loop over input quantitative maps
    for j = 1:nv
        vals = vals(ind + (j-1)*np);
        nval = numel(vals);
        Mstr = M_label{j,1};

        % Loop over functions
        fcn_ind = M_label{j,2};
        for k = 1:numel(fcn_ind)
            switch fcn_ind(k)
                case 1 % Mean
                    res = mean(vals);
                    res_label = {[Mstr,"_mean"]};
                case 2 % Percent volume
                    res = nnz(vals)/nval *100;
                    res_label = {[Mstr,"_%"]};
                case 3 % Percent volume for each value (PRM)
                    nprm = numel(prm_val{j});
                    res = zeros(1,nprm);
                    for i_prm = 1:nprm
                        res(i_prm) = (vals==prm_val{iprm})/nval *100;
                    end
                    res_label = Mstr + "_" + prm_val{j};
            end
            Vind = ii + (1:numel(res));
            V(i,Vind) = res;
            V_label(Vind) = res_label;
        end

    end

end

%% Replace zero value with average between adjacent values











uStrahl = unique(B(:,6));
for i = 2:numel(uStrahl)
    ind = find(B(:,6)==uStrahl(i));

end





function P = CTPerf(exp,ins,J,mask,voxsz)

% Initialize parameters
ksp = 30; % mm
se = strel('cuboid',[7,7,3]);
d = size(exp);
fov = d.*voxsz;

% Naive perf 
P = abs(exp - ins.*J);

% Determine subregion centers:
[gx,gy,gz] = meshgrid(round((ksp/2:ksp:fov(1))/voxsz(1)),...
                      round((ksp/2:ksp:fov(2))/voxsz(2)),...
                      round((ksp/2:ksp:fov(3))/voxsz(3)));
K = length(gx);

for sub_i = 1:K(1)
    ind_i = (sep_i(sub_i)+1) : sep_i(sub_i+1);
    for sub_j = 1:K(2)
        ind_j = (sep_j(sub_j)+1) : sep_j(sub_j+1);
        for sub_k = 1:K(3)
            ind_k = (sep_k(sub_k)+1) : sep_k(sub_k+1);
            
            % Isolate subregion:
            t_exp = exp(ind_i,ind_j,ind_k);
            t_ins = ins(ind_i,ind_j,ind_k);
            t_mask = mask(ind_i,ind_j,ind_k);
            
            t_P = CTperf_sub(t_exp,t_ins,t_mask,voxsz);
    
            tau = 
        end
    end
end

end

function P = CTperf_sub(exp,ins,mask,voxsz)
    % Find voxel locations:
    N = nnz(tmask);
    x = ind2sub(size(tmask),find(tmask)) * diag(voxsz) - repmat(voxsz/2,N,1);

    % Find knot locations:
    % Paper uses dart throwing: 
    %       Dart Throwing on Surfaces, Cline D, et al, Comp. Graph. Forum, 2009
    z = [];
    L = length(z);

    % Calculate C_ij
    sig = 1;
    C = zeros(N,L);
    for j = 1:L
        C(:,j) = sum((x-repmat(z(j,:),N,1).^2),2);
    end
    C = exp(-sig*C);
    C = C./repmat(sum(C,2),1,L);

    Ahat = A * C;


    p = lsqlin(C/
end
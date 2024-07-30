function p = tree_prep_pipe(pID,fn_seg,fn_airways,outD)
p = struct('pixDim',{[]},'lobe_ids',{''},'B',{[]},'N',{[]},'lobe_surfs',{[]});
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
%% 0. load data (ariways and lobes)

% Lobe Segmentation
L = double(niftiread(fn_seg));
dim = size(L);

% Airways
A_info = niftiinfo(fn_airways);
A = logical(niftiread(A_info));
p.pixDim = A_info.PixelDimensions;

% Check that marix size matches
if any(size(A)~=dim)
    error('Matrix dimensions do not match.');
end

% Check lobe IDs
lobe = getLobeTags(L);
lobe = lobe(cellfun(@(x)numel(x)==1,{lobe.val})); % Use singular regions

%% 1. the major airway centreline graph

[A_D,A_D_idx] = bwdist(~A,'euclidean'); % to use for radii approximation

A_vox = im2coords3Dbin(A);
A_cl = bwskel(A);
A_cl_vox = get_M(A_cl,A_D,A_D);
A_cl_vox = A_cl_vox(:,1:4);

hf = figure('Name','Centreline Voxels'); ha = axes(hf);
plot3(ha,A_cl_vox(:,1),A_cl_vox(:,2),A_cl_vox(:,3),'.');
axis(ha,'equal');
view(ha,[0,1,0]);
saveas(hf,fullfile(outD,[pID,'.cl_vox.fig']));
close(hf)

[~,node,link] = Skel2Graph3D(A_cl,0);

hf = figure('Name','Centreline Graph'); ha = axes(hf);
plot_graph(ha,node,link)
A_bound = boundary(A_vox,1); % N.B. not really tight enough, but low priority
A_tri = triangulation(A_bound,A_vox);
trisurf(A_tri,'FaceColor','black','FaceAlpha',0.1,'EdgeColor','none');
grid(ha,'on');
axis(ha,'equal');
view(ha,[0 1 0]);
saveas(hf,fullfile(outD,[pID,'.cl_graph.fig']));
close(hf);


%% 2. graph link radii assignment
%
n = size(link,2); R = zeros(1,n);
d2 = @(v1,v2) sqrt(sum((v1-v2).^2));
for i = 1:n
    brc = link(i);
    pnts = brc.point;
    len = length(pnts);
    
    c = (len+1)/2;
    lim = round([c-len/4 c+len/4]);
    pnts = pnts(lim(1):lim(2));
    
    tmp = zeros(1,length(pnts));
    
    for j = 1:length(tmp)
        [x,y,z] = ind2sub(dim,pnts(j));
        CLpnt = [p.pixDim(1)*x p.pixDim(2)*y p.pixDim(3)*z]; % point on centreline
        ind = A_D_idx(x,y,z);
        [x,y,z] = ind2sub(dim,ind);
        NRpnt = [p.pixDim(1)*x p.pixDim(2)*y p.pixDim(3)*z]; % nearest point from bwdist
        tmp(j) = d2(CLpnt,NRpnt);                   % Euclidean distance in mm between points
    end
    
    R(i) = mean(tmp); % approximate radius in mm from mean of sampled midpoints
end


%% 3. formatting data into branches and nodes, more familiar structures
%
s = size(link,2);
B = zeros(s,3); % node 1, node 2, radius
for i = 1:s
    B(i,1) = link(1,i).n1;
    B(i,2) = link(1,i).n2;
    B(i,3) = R(i);
end
% NOD_ids=unique([B(:,1);B(:,2)]);
% NOD=zeros(length(NOD_ids),4);
% for i=1:length(NOD_ids)
%     N_id=NOD_ids(i);
%     NOD(i,1)=N_id;
%     NOD(i,2)=node(1,N_id).comx;
%     NOD(i,3)=node(1,N_id).comy;
%     NOD(i,4)=node(1,N_id).comz;
% end
N_ids = unique([B(:,1);B(:,2)]);
N = zeros(max(N_ids),4);
for i = 1:numel(N_ids)
    N_id = N_ids(i);
    N(N_id,1)=N_id;
    N(N_id,2) = node(1,N_id).comx;
    N(N_id,3) = node(1,N_id).comy;
    N(N_id,4) = node(1,N_id).comz;
end


%% 4. tree orientation

% 1. distal nodes
n = size(N,1);
distN = zeros(n,1);
for i = 1:n
    chk = B(:,1:2)==N(i,1);
    if sum(sum(chk)) == 1
        distN(i) = 1;
    end
end
distN = N(logical(distN),:);

% 2. distal branches
n = size(B,1);
distB = zeros(n,1);
for i = 1:n
    if ismember(B(i,1),distN(:,1)) || ismember(B(i,2),distN(:,1))
        distB(i) = 1;
    end
end
B2 = [(1:n)' B];
distB = B2(logical(distB),:);
maxR = max(distB(:,end)); % largest radius
trach = B2(B2(:,4)==maxR,:); % root branch (assuming trachea)
idT = trach(1); % row index for trachea
B = [B(idT,:) ; B(1:idT-1,:) ; B(idT+1:end,:)]; % place root first in B

% 3. check root is oriented proximal to distal
if sum(sum(B(:,1:2)==B(1,2)))==1 % check if external node is 'distal' here
    B(1,:) = [B(1,2) B(1,1) B(1,3)]; % transpose if required
end

% 4. orient tree relative to root branch
rootN = B(1,1); % root node
G = graph(B(:,1),B(:,2));
TR = shortestpathtree(G,rootN); e=TR.Edges;
B2 = e.EndNodes;

hf = figure; ha = axes(hf);
hold(ha,"on");
for i = 1:size(B2,1)
    try
        pnt = [N(B2(i,1),2) N(B2(i,1),3) N(B2(i,1),4)];
        vec = [N(B2(i,2),2)-N(B2(i,1),2) N(B2(i,2),3)-N(B2(i,1),3) ...
               N(B2(i,2),4)-N(B2(i,1),4)];
        quiver3(ha,pnt(1),pnt(2),pnt(3),vec(1),vec(2),vec(3),'b-',...
            'LineWidth',2,'MaxHeadSize',0.5);
    catch err
        disp(getReport(err))
        continue
    end
end
saveas(hf,fullfile(outD,[pID,'.oriented_tree.fig']));
close(hf);

B2 = B2(end:-1:1,:); % place root at top, tree is now oriented

% 5. include radii from B in oriented tree
for i = 1:size(B2,1)
    brc = B2(i,1:2); % prox id dist id for branch
    [a,b] = ismember(brc,B(:,1:2),'rows');
    if a
        B2(i,3) = B(b,3); % bring over radii
    else
        [~,b] = ismember([brc(2) brc(1)],B(:,1:2),'rows');
        B2(i,3) = B(b,3);
    end
end
B=B2; % complete conversion, B is now an oriented branch matrix with NOD as
% associated matrix of nodes, with the root (trachea) as first element


%% tree generation, Strahler and Horsfield order assignment
%

H=calc_ho(B); H=H(:,end); S=calc_so(B); S=S(:,end);

% calculating generation no. given a recursive function defined in footer
n = size(B,1);
G = zeros(n,1);
for i = 1:n
    G(i) = gen_no(B(i,:),B);
end
B = [B G S H]; % append generation, Strahler and Horsfield (in that order)

hf = figure; t = tiledlayout(hf,1,3);
for i = 4:6
    plot_tree(nexttile(t),B,N,i);
end
saveas(hf,fullfile(outD,[pID,'.raw_tree_orders.fig']));
close(hf);

hf = figure; t = tiledlayout(hf,1,3);
T={'Generation' 'Strahler' 'Horsfield'};
for i=4:6
    ha = nexttile(t);
    h = histogram(ha,B(:,i));
    h.BinWidth = 1;
    title(ha,T{i-3})
end
saveas(hf,fullfile(outD,[pID,'.raw_tree_hists.fig']));
close(hf);

%% convert nodes to mm and calculate branch lengths
%
% 0. add branch lengths to end (7) of B
%
% first converting NOD to mm
for i = 2:4
    N(:,i) = N(:,i)*p.pixDim(i-1);
end

% now branch length assignment using d2 metric
n=size(B,1); Lengs=zeros(n,1);
for i = 1:n
    n1 = N(B(i,1),:);
    n2 = N(B(i,2),:);
    Lengs(i) = d2(n1(2:end),n2(2:end));
end
B = [B Lengs];

hf = figure; t = tiledlayout(hf,1,2);
ha = nexttile(t);
histogram(ha,B(:,3));
title(ha,'Radii')
xlabel(ha,'airway radius (mm)')
ha = nexttile(t);
histogram(ha,B(:,end));
title(ha,'Lengths')
xlabel(ha,'airway length (mm)')
saveas(hf,fullfile(outD,[pID,'.pre_QC_lengths_and_radii.fig']));
close(hf);

%% branch QC algorithm
%

% Idea (Brody): using thresholds t1 and t2 (set below), remove branches
% through comparison of their length (L) to parent branch length (Lp)

% L < t1*Lp => remove
% t1*Lp <= L <= t2*Lp => keep (+ match length with sibling <test>)
% L > t2*Lp => accept

% After each pass, a clean up step will follow removing interm. nodes
% Process will conintue until no terminal branches flagged for alteration

B2=B; t1=0.1; t2=0.4; flg=1;
L_fix = 2; % addition to algorithm, a fixed minimum length (mm) for all branches

while flg % whole process loop
   
    n = size(B2,1); elim = ones(n,1); flg = 0;
    
    % step 1: testing terminal branches to alter based on length
    for i = 2:n % always skip root
        b = B2(i,:);
        if not(ismember(b(2),B2(:,1))) % test if terminal
            B_P = B2(B2(:,2)==b(1),:); % obtain parent branch
            L_P = B_P(end); % parent branch length
            if b(end)<t1*L_P || b(end)<L_fix
                flg = 1; % put flag back up to continue trimming
                elim(i) = 0; % mark for elimination
            elseif not(b(end)>t2*L_P)
                % take action on assumption of 'truncated branch'
            end
        end
    end
    
    B2 = B2(logical(elim),:); % remove marked branches
    
    % cleaning of intermediate nodes...
    % idea check all branches except the root for number of children, if
    % that number is 1, merge parent and child into one new branch,
    % obtaining new properties as averages of radii and addition of lengths
    
    flg_clean = 1;
    while flg_clean 
        n=size(B2,1); elim=ones(n,1); flg_clean=0;
        for i = 2:n
            b = B2(i,:);
            ch = B2(B2(:,1)==b(2),:);
            if size(ch,1)==1 % branch has one child
                [~,j] = ismember(ch,B2,'rows');
                % i and j now index rows of B2 to replace with merged brch
                B2(i,:) = [b(1) ch(2) 0.5*(b(3)+ch(3)) b(4:6) (b(7)+ch(7))];
                elim(j) = 0; 
                B2 = B2(logical(elim),:);
                flg_clean = 1;
                break
            end
        end
    end
end

%
B = B2;
% re-evaluate Strahler and Horsfield orders
H=calc_ho(B); H=H(:,end); S=calc_so(B); S=S(:,end);
B(:,5) = S; B(:,6) = H;
n = size(B,1);
G = zeros(n,1);
for i = 1:n
    G(i) = gen_no(B(i,:),B);
end

B(:,4) = G; % append generation, Strahler and Horsfield (in that order)
% remove redundant nodes from NOD?
% Leaving in for now to maintain row ID agreement


%% obtaining lobe surface maps in mm (to match tree)
%
nl = numel(lobe);
L_surfs = cell(1,nl);
for i = 1:nl
    tmp = L==lobe(i).val; % isolate lobe
    tmp = im2coords3Dbin(tmp); % convert to [x,y,z]
    
    for j = 1:3
        tmp(:,j) = p.pixDim(j)*tmp(:,j);
    end
    
    s = size(tmp);
    tmp2 = datasample(tmp,floor(s(1)/10)); % random sample of 10%
    tmp3 = boundary(tmp2); % obtain indices for boundary set of voxels
    TR = triangulation(tmp3,double(tmp2));
    L_surfs{i} = TR;
end


%% final visual check
%
hf = figure('Name','Final Visual'); ha = axes(hf);
plot_tree(ha,B,N,4)
for i=1:5
    trisurf(L_surfs{i},'FaceColor','b','FaceAlpha',0.05,'EdgeColor','none');
end
saveas(hf,fullfile(outD,[pID,'.final_visual.fig']));
close(hf);


%% saving data to pass to Brody / use for airway growing algorithm
%
% basic approximation given by multiplying by average of pixel dims, n.b.
% typically z dim is different to x and y dims which themselves match

p.lobe_ids = [lobe.val];
p.lobe_surfs = L_surfs;
p.B = B;
p.N = N;
save(fullfile(outD,[pID,'.tree_data.mat']),'-struct','p');

warning('on','MATLAB:triangulation:PtsNotInTriWarnId')

% generation no., given all roads lead to root (B(1,:)).
function y = gen_no(b,B) % given branch b from B, output generation no.
    if b(1,1)==B(1,1) % root, define as generation 1
        y=1;
    else % recursion using parent branch
        b=B(B(:,2)==b(1,1),:);
        y=1+gen_no(b,B);
    end
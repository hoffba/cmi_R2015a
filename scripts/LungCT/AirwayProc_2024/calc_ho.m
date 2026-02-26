function HO = calc_ho(B)
% calculate Horsfield order from branches B

% Source: https://en.wikipedia.org/wiki/Strahler_number
% If the node is a leaf (has no children), its HO number is one.
% Otherwise its ho number is 1 greater than the largest ho of its children.
% N.B. copied from SO script, almost identical

nB = size(B,1);
HO = zeros(nB,1);
Nprox = B(:,1);
Ndist = B(:,2);

% Find terminal branches
% - any branch whose distal node does not occur in the list of proximal nodes
activeB = findTerminalBranches(B);
HO(activeB) = 1;

% Start loop at terminal branches' parents
activeB = find(ismember(Ndist,Nprox(activeB)));

% Assign HO values starting at terminal branches until all branches are assigned
stat = true;
while any(stat) && ~isempty(activeB)

    % Try to assign value to active branches
    naB = numel(activeB);
    stat = false(naB,1);
    for i = 1:naB
        % Find values of children
        val = HO(Nprox==Ndist(activeB(i)));
        if all(val)
            HO(activeB(i)) = max(val)+1;
            stat(i) = true;
        end
    end

    % Add parents of branches with set HO values
    % Remove branches with set HO
    addB = find(ismember(Ndist,Nprox(activeB(stat))));
    activeB = [activeB(~stat);addB];
end

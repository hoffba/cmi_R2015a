function SO = calc_so(B)
% calculate Strahler order from branch marix B2 generated from calcNOD

% Source: https://en.wikipedia.org/wiki/Strahler_number
% If the node is a leaf (has no children), its Strahler number is one.
% If the node has one child with Strahler number i, and all other children have Strahler numbers less than i, then the Strahler number of the node is i again.
% If the node has two or more children with Strahler number i, and no children with greater number, then the Strahler number of the node is i + 1.

nB = size(B,1);
SO = zeros(nB,1);
Nprox = B(:,1);
Ndist = B(:,2);

% Assign terminal branches a value of 1
activeB = ~ismember(Ndist,Nprox);
SO(activeB) = 1;

% Start searching parents of terminal branches
activeB = find(ismember(Ndist,Nprox(activeB)));

stat = true;
while any(stat) && ~isempty(activeB)

    % Try to assign value to active branches
    naB = numel(activeB);
    stat = false(naB,1);
    for i = 1:naB
        % Find values of children
        val = SO(Nprox==Ndist(activeB(i)));
        if all(val)
            uval = unique(val);
            if isscalar(uval)
                val = uval + 1; % All children have the same SO
            else
                val = max(val);
            end
            SO(activeB(i)) = val;
            stat(i) = true;
        end
    end

    % Add parents of branches with set HO values
    % Remove branches with set HO
    addB = find(ismember(Ndist,Nprox(activeB(stat))));
    activeB = [activeB(~stat);addB];
end
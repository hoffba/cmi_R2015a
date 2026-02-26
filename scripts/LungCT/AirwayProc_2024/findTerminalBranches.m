function iBterm = findTerminalBranches(B)
% Finds terminal branches of tree
% Assumes branch nodes are ordered proximal (B(:,1)) to distal (B(:,2))
% Returns indices of terminal branches in the branch list

iBterm = find(~ismember(B(:,2),B(:,1)));

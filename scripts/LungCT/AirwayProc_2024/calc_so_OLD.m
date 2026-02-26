function y = calc_so(B)
% calculate Strahler order from branch marix B2 generated from calcNOD

% Source: https://en.wikipedia.org/wiki/Strahler_number
% If the node is a leaf (has no children), its Strahler number is one.
% If the node has one child with Strahler number i, and all other children have Strahler numbers less than i, then the Strahler number of the node is i again.
% If the node has two or more children with Strahler number i, and no children with greater number, then the Strahler number of the node is i + 1.

s=size(B); n=s(1); % no. of branches
SO=zeros(n,1); % column of SO numbers
B2=[B (1:n)']; % append ID column

% SO1 assignment, criteria: any branch who's distal node does not occur in
% the list of proximal nodes

for i=1:n
    if ~(ismember(B2(i,2),B2(:,1)))
        SO(i)=1; % terminal branch
    end
end

% now shall repeatedly sweep tree assigning SO until all SO is assigned
while ismember(0,SO)
    for i=1:n
        if SO(i)==0 % branch still needs an SO assigned
            % get list of all daughter branches
            tmp=B2(B2(:,1)==B2(i,2),:);
            s=size(tmp); flg=0; %flag is a daughter is unassigned
            SOchk=zeros(1,s(1)); %record SOs of daughters
            for j=1:s(1)
                if SO(tmp(j,end))==0
                    flg=1; %daughter yet to be assigned SO
                    break
                end
                SOchk(j)=SO(tmp(j,end));
            end
            if ~flg % all daughters have an SO, now to assign new SO
                if sum(ismember(SOchk,max(SOchk)))>1
                    SO(i)=max(SOchk)+1;
                else
                    SO(i)=max(SOchk);
                end
            end
        end
    end
end
y=[B2(:,1:3) SO]; % provide B3 with key information
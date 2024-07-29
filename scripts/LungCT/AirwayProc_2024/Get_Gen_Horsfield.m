function [Gen, Horsfield, Strahler] = Get_Gen_Horsfield(Branches, lower_branches)
%% Function to calculate Generation, Horsfield, and Strahler for all branches
% in a network

Horsfield = zeros(length(Branches),1);
Gen = Horsfield; 
Strahler = Gen; 

%% Calculate Horsfield and Strahler order
parents = NaN(length(Branches),1);
daughters = cell(length(Branches),1);
for i = 2:length(Branches)
    parents(i) = find(Branches(:,3) == Branches(i,2));
    daughters{i} = find(Branches(:,2) == Branches(i,2));
end


Horsfield(lower_branches) = 1;
Strahler(lower_branches) = 1; 
for i = 1:length(lower_branches)
    current_branch = lower_branches(i);
   
    while current_branch ~= 1
        new_parent = parents(current_branch);
        new_daughters = daughters{current_branch};

        maxH = max(Horsfield(new_daughters));
        Horsfield(new_parent) = maxH + 1; 
        curS = Strahler(new_daughters); 
        maxS = max(curS); 

        if sum(curS == maxS) > 1
            Strahler(new_parent) = maxS + 1; 
        else
            Strahler(new_parent) = maxS; 
        end        
        current_branch = new_parent;
    end
end

%% Calculate generation (recursively)
Gen(1) = 1;

current_branch = 1; 
Gen = Calculate_Gen(Branches, Gen, current_branch, lower_branches); 

end

function Gen = Calculate_Gen(Branches, Gen, current_branch, lower_branches)

branches_below = find(Branches(:,2)==Branches(current_branch,3));

current_gen = Gen(current_branch);

for i = 1:length(branches_below)
   Gen(branches_below(i)) = current_gen + 1;
   
   Gen = Calculate_Gen(Branches, Gen, branches_below(i), lower_branches);
end
end
        
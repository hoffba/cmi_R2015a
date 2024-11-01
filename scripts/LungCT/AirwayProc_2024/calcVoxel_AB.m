% Function to assign ct-based measurement to distal point in different ways

% 09.01.2022
% Inputs...
% brc ~ [proximal node, distal node, radius, lobe
% bro ~ [x,y,z] ~ nodes
% imb ~ imbio n x 4 matrix, n ~ order 10^6or7 4th col. ~ functional val.
% so ~ strahler order for branches (indexed same as brc)
% e.g. delHU
% px ~ pixel dimension

% Outputs...
% brc: expHU and InsHU for all nodes

function out = calcVoxel_AB(brc,bro,imb,so,r,px,b,Nod)

%% 1. Find terminal branches(sorting by soc)

n = length(so); % no. of branches to go through

soCnt = length(unique(so)); % number of SO, assuming 1,2,...,soCnt
soCnts = zeros(1,soCnt);

for i = 1 : soCnt
   
    soCnts(i) = sum(so==i); % number of each SO occuring
    
end

brc = [brc so]; % attach so info to 5th column of brc
brc = sortrows(brc,5); % arrange with lowest SO at top
brc = [brc nan(n,17)]; % fill several col (to record other measurements) with NaN

%% 2. Creating 3D matrix to store Mask and real CT data
imData = nan(max(imb(:,1)),max(imb(:,2)),max(imb(:,3))); % create 3D matrix to store the raw CT data and mask  
eHUData = imData; % to store expHU
iHUData = imData; % to store insHU

broZ = round(bro./px); %% Convert the coordinate into intergral
for i = 1:size(imb,1)
    eHUData(imb(i,1),imb(i,2),imb(i,3)) = imb(i,4); % Store the data 
    iHUData(imb(i,1),imb(i,2),imb(i,3)) = imb(i,5);
end

se = calcSe(px,r); % n * 3 matrix, the distance between the nodes within the range and the central nodes(end of the terminal branch).

%% 3. Voxel Association

len = size(se,1);
for i = 1:soCnts(1)
    index = se + [broZ(brc(i,2),1),broZ(brc(i,2),2),broZ(brc(i,2),3)]; % get the coordinate of all possible nodes (no mather if there have some values)
    num = 1; 
%     tmpIndex = zeros(len,3); % Used to store the coordinate of nodes within the range R.
    tmpData = NaN(len,2); %[eHU iHU delHU]
    for j = 1:len
        if index(j,1) > 0 && index(j,2) >0 && index(j,3) >0 
            if index(j,1)<=size(eHUData,1) && index(j,2)<=size(eHUData,2) && index(j,3)<=size(eHUData,3)  % make sure it stays in the boundary 
                if ~(isnan(eHUData(index(j,1),index(j,2),index(j,3))))     % whether there is data at the point
                    tmpData(num,1) = eHUData(index(j,1),index(j,2),index(j,3));
                    tmpData(num,2) = iHUData(index(j,1),index(j,2),index(j,3));              
%                     tmpIndex(num,:) = index(j,:);
                    num = num+1;                
                end
            end
        end  
    end
    tmpData = tmpData(1:(num-1),:);

%     cellIndex{i} = tmpIndex(1:(num-1),:);   % uncommand these lines if
%     wish to return the coordinate of nodes with in the range
%     cellIndex{i} = cat(2,tmpIndex(1:(num-1),:),tmpData); %  

    if ~isnan(tmpData)
        brc(i,6) = 0; %% indicate if it is real branches; 0 means 'no', 1 means 'yes'
        brc(i,7) = mean(tmpData(:,1),1); %eHU
        brc(i,8) = mean(tmpData(:,2),1); %iHU
        brc(i,9) = sum(tmpData(:,1) > -856); %% Numbers of nodes beyond -856
        brc(i,10)= num - 1 - brc(i,9); %% Numbers of nodes below -856
        brc(i,11)= sum(tmpData(:,2) > -950); %% Numbers of node beyond -950
        brc(i,12)= num - 1 - brc(i,11); %% Numbers of node below -950
        
        
        %%% PRM Classification
        brc(i,13) = sum(tmpData(:,2) > -810); %% p_a
        brc(i,14) = sum(tmpData(:,1) < -856 & tmpData(:,2) > -950 & tmpData(:,2) < -810); %% y_a/ fsad
        brc(i,15) = sum(tmpData(:,1) > -856 & tmpData(:,2) > -950 & tmpData(:,2) < -810); %% g_a/ 
        brc(i,16) = sum(tmpData(:,1) < -856 & tmpData(:,2) < 950); %% r_a / prm_emph
        
        %%% Result of Voting 
        brc(i,17) = sum(tmpData(:,1) < -856) / (num-1); % percentage of nodes that are gas trapping
        brc(i,18) = sum(tmpData(:,2) < -950) / (num-1); % percentage of nodes that are emphsema
        brc(i,19) = brc(i,13) / (num-1); %% prm classification
        brc(i,20) = brc(i,14) / (num-1);
        brc(i,21) = brc(i,15) / (num-1);
        brc(i,22) = brc(i,16) / (num-1);
    
    else
        brc(i,6:22) = zeros(1,17); %% Boxes with no voxels are signed with zero
    end
end

%% 4. Replace those nodes assigned to zero with mean value
indZ = brc(:,7)==0;
indNz = brc(1:soCnts(1),7)~=0;
brc(indZ,7:22) = mean(brc(indNz,7:22), 1) .* ones(1,nnz(indZ))';

%% 5. Calculate the value around other branches with higer so
len = soCnts(1);
brc = sortrows(brc,5);
for i = 2:soCnt  %% sort each row so as to make it easier to find mother branches
    brc((len+1):(len+soCnts(i)),:) = sortrows(brc((len+1):(len+soCnts(i)),:),2,'descend');
    len = len + soCnts(i);   
end


%% 6. Indentify whether the branches belong to real tree
trn = NaN(size(Nod,1),2); %% [row_id in B,   row_id in BRC]
for rows = 1:size(Nod,1)
    trn(rows,1) = rows;
    ind_tmp = find(bro(:,1)==Nod(rows,2) & bro(:,2)==Nod(rows,3)& bro(:,3)==Nod(rows,4));
    if ~isnan(ind_tmp) % if the nodes exist in the conducting airway
        trn(rows,2) = ind_tmp;
    end
end
sim_id = [trn(b(:,1),2) trn(b(:,2),2)]; % row_id in BRC for all the real branches

%% 7. Calc all the nodes with Strahler order larger than 1 and record the position of missing value
k = 0;
reChk = [1 0 0 0 0 0 0] .*ones(200,1);  %% Store the branches whose daughter branches haven't been assign with value yet.
% [postID, distID, rowIndMotherNode, rowInd_DaughterNode1, rowInd_DaughterNode2, rowInd_DaughterNode3]
    for i = soCnts(1)+1:size(brc,1)
        if ~isnan(find(sim_id(:,1) == brc(i,1) & sim_id(:,2) == brc(i,2))) % If it is real airway
            brc(i,6) = 1;
        else
            brc(i,6) = 0;
        end
        tmpi = find(brc(:,1)==brc(i,2)); %% row_th of daughter nodes
        tmpa = brc(tmpi,7:end); %% values of daughter nodes
        %%% Finding lobe_id
        if var(brc(tmpi,4)) == 0 % Lobe id of all daughter nodes are equal
            brc(i,4) = brc(tmpi(1),4);
        end
        %%% Assign value
        if sum(isnan(tmpa)) == 0 %~isnan(tmpa(1)) && ~isnan(tmpa(2)) which means both two daughter nodes have been found
            Nn = (tmpa~=0); 
            if sum(Nn(:,1)) == 0 % if daughter nodes don't have value
                brc(i,7:end) = zeros(1,10);
            else
                brc(i,7:8) = mean(tmpa((Nn(:,1)~=0),1:2),1);
                %%% 7th: expHU, 8th: insHU
                brc(i,9:16)= sum( tmpa(Nn(:,1)~=0, 3:10) ,1); 
                %%% 9-10th: Nums of Gt/Non_GT nodes, 11-12th: Nums of Emph/Non-Emph nodes 
                %%% 13-16th: Nums of PRM Classification
                brc(i,17:22)= mean( tmpa(Nn(:,1)~=0, 11:16) ,1); 
                %%% 17-18th: Mean of Pct_of_Gt/Pct_of_Emph in daugther nodes
                %%% 19-22th: Mean of Pct_PRM_Classification
            end
        else 
            reChk(k+1,1:3) = [brc(i,1:2) i]; %% [prox.id, dist.id, rows] 
            reChk(k+1,4:3+size(tmpi,1)) = tmpi'; %% dist.daughers'row (2to3 daughters),
            k = k+1;
        end
    end

%% 8. Recheck those missing data by changing soring order of col3 each time.
reChk = reChk(1:k,:);
a = 1;
dir = {'descend','ascend'}; k = 1; % change the direction of finding nodes
%%% reChk = [prox.id, dist.id, row_id_Mother, row_id_daughter1, row_id_daughter2, row_id_daughter3, chk]
while a >0 %% whether there is NaN in reChk
    ii = find(reChk(:,6)>0); % find which nodes have three daughter branches
    for i = 1:size(reChk,1)
        if ismember(i,ii) % Indicate a node with three daughter nodes
            tmpi = reChk(i,4:6); % tmpi are the row-th for daughter nodes
        else
            tmpi = reChk(i,4:5);
        end
        if sum(isnan(brc(tmpi,7:22))) == 0 %% if tmpi has some value
            if var(brc(tmpi,4)) == 0
                brc(reChk(i,3),4) = brc(tmpi(1),4);
            end
            brc(reChk(i,3),7:8) = mean(brc(tmpi,7:8),1);
            brc(reChk(i,3),9:16) = sum(brc(tmpi,9:16),1);
            brc(reChk(i,3),17:22) = mean(brc(tmpi,17:22),1);
            reChk(i,7) = 1;
        end
        
    end
    a = sum(isnan(brc(:,7:22)));
    k = mod(k,2)+1;
    reChk = sortrows(reChk,3,dir(k));     
end






out = brc;
% out{2} = cellIndex;

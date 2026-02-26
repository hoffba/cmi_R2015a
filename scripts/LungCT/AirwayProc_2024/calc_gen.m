function [G,B] = calc_gen(B,activeN)    
% Calculate branch generations

if nargin<2
    activeN = B(1,1);
end

nB = size(B,1);
G = zeros(nB,1);
gi = 0;
while any(~G)

    % Find branches
    Bi = find(any(ismember(B(:,1:2),activeN),2) & ~G);

    % Increment generation and set values
    gi = gi+1;
    G(Bi) = gi;
    
    % Loop over branches to correct direction
    for i = 1:numel(Bi)
        % Check if prox/dist need switch
        if ~ismember(B(Bi(i),1),activeN)
            B(Bi(i),1:2) = B(Bi(i),[2,1]);
        end
    end

    % update active node list
    activeN = B(Bi,2);

end

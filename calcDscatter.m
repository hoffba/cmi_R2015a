function [D,Dmax] = calcDscatter(X,Y)
% calcDscatter returns the population density index vector for the 2D X-Y plot
%  -- simplified from DSCATTER

% Reference:
% Paul H. C. Eilers and Jelle J. Goeman
% Enhancing scatterplots with smoothed densities
% Bioinformatics, Mar 2004; 20: 623 - 628.

minx = min(X,[],1);
maxx = max(X,[],1);
miny = min(Y,[],1);
maxy = max(Y,[],1);

nbins = [min(numel(unique(X)),200) ,min(numel(unique(Y)),200) ];
lambda = 20;

edges1 = linspace(minx, maxx, nbins(1)+1);
edges1 = [-Inf edges1(2:end-1) Inf];
edges2 = linspace(miny, maxy, nbins(2)+1);
edges2 = [-Inf edges2(2:end-1) Inf];

[n,~] = size(X);
bin = zeros(n,2);
% Reverse the columns to put the first column of X along the horizontal
% axis, the second along the vertical.
[~,bin(:,2)] = histc(X,edges1);
[~,bin(:,1)] = histc(Y,edges2);
H = accumarray(bin,1,nbins([2 1])) ./ n;
G = smooth1D(H,nbins(2)/lambda);
F = smooth1D(G',nbins(1)/lambda)';

Dmax = max(F(:));
F = F./Dmax;
ind = sub2ind(size(F),bin(:,1),bin(:,2));
D = F(ind);

function Z = smooth1D(Y,lambda)
[m,~] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;


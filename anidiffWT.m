function img = anidiffWT(img,wttype,lvls,niter,dt,kappa,opt)
% Performs 2D image filtering using anisotropic diffusion in wavelet domain
% Syntax:
%   img = anidiff(img,);
% Inputs:

% Wavelete transform:
[C,S] = wavedec2(img,lvls,wttype);

% Anisotropic diffusion filtering of detail images:
tS = prod(S,2);
N = tS.*[1;3*ones(lvls,1);0];
for il = 1:lvls
    for ii = 1:3
        
        % Reshape data to 2D:
        ind = sum(N(1:il)) + (ii-1)*tS(il+1) + (1:tS(il+1));
        tC = reshape(C(ind),S(il+1,:));
        
        % Filter the detail image:
        tC = anisodiff2D(tC,niter,dt,kappa,opt);
        
        % Replace original detail data:
        C(ind) = tC;
        
    end
end

% Wavelet reconstruction:
img = waverec2(C,S,wttype);
    
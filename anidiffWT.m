function img = anidiffWT(img,wttype,lvls,niter,dt,kappa,opt)
% Performs 2D image filtering using anisotropic diffusion in wavelet domain
% Syntax:
%   img = anidiff(img,);
% Inputs:

% Wavelete transform:
[C,S] = wavedec2(img,lvls,wttype);

% Anisotropic diffusion filtering:
tS = prod(S,2);
N = tS.*[1;3*ones(lvls,1);0];
for il = 1:lvls
    for ii = 1:3
        ind = sum(N(1:il)) + (1:tS(il+1));
        C(ind) = anisodiff2D(reshape(C(ind),S(il+1,:)),niter,dt,kappa/(lvls-il+1),opt);
    end
end

% Wavelet reconstruction:
img = waverec2(C,S,wttype);
    

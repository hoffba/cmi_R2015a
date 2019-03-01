function batch_fluxLungReg(username,catalog_name)
% Input:
%   catalog_name = name of text file cataloging lung data to process

% Set Elastix options
elxbindir = '/nfs/turbo/umms-cgalban/Elastix';
elxfname = fullefile(elxbindir,{'ElastixParameters_0_affine.txt','ElastixParameters_1_warp.txt'});
opts = {'f_clim',[-inf,0],...
        'm_clim',[-inf,0],...
        'f_dilateN',5*ones(1,3),...
        'f_filtN',ones(1,3),...
        'm_filtN',ones(1,3),...
        'jac',true};
    
% Read catalog file:    

% Loop over cases:
for i = 1:N
    odir = '';
    [~,workingdir] = fileparts(odir);
    % Temporarily use scratch storage for speed:
    workingdir = fullfile('/scratch/cgalban_fluxod',username,workingdir);
    f_name = '';
    fv_name = '';
    m_name = '';
    mv_name = '';
    batch(@autoRegLung,0,[{elxbindir,odir,workingdir,elxfname,f_name,fv_name,m_name,mv_name,},opts]);
end
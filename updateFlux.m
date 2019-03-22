function updateFlux(username)

fList = matlab.codetools.requiredFilesAndProducts('batch_fluxLungReg.m')';
system(strcat('scp',sprintf(' %s',fList{:}),...
    sprintf(' %1$s@flux-xfer.arc-ts.umich.edu:/home/%1$s/Matlab',username)));
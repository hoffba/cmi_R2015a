function J = def2jac(D,voxsz)
J = [];
d = size(D);
if (length(d)==4) && (d(4)==3)
    J = zeros([d(1:3),9]);
%     D = padarray(D,ones(1,3),'replicate');
    for i = 1:3
        fprintf('Calculating gradient %u of 3 ... ',i);
        j = (i-1)*3;
        [FY,FX,FZ] = gradient(D(:,:,:,i),voxsz(1),voxsz(2),voxsz(3));
        J(:,:,:,j+1) = FX;
        J(:,:,:,j+2) = FY;
        J(:,:,:,j+3) = FZ;
        J(:,:,:,j+i) = J(:,:,:,j+i)+1;
        fprintf('done\n');
    end
end
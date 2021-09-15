function [L1,L2,L3] = EigenMatrix3x3M(T)
%% Matlab version

image_size = size(T(:,:,:,1));  
 
x = [reshape(T(:,:,:,1), 1,1,[]), reshape(T(:,:,:,2), 1,1,[]), reshape(T(:,:,:,3), 1,1,[]);...
    reshape(T(:,:,:,4), 1,1,[]), reshape(T(:,:,:,5), 1,1,[]), reshape(T(:,:,:,6), 1,1,[]);...
    reshape(T(:,:,:,7), 1,1,[]), reshape(T(:,:,:,8), 1,1,[]), reshape(T(:,:,:,9), 1,1,[])];
L = eig3(x);
L1 = reshape(L(1,:), image_size);
L2 = reshape(L(2,:), image_size);
L3 = reshape(L(3,:), image_size);
% L1 = zeros(size(x1));
% L2 = zeros(size(x1));
% L3 = zeros(size(x1));
% s2 = size(x1,2);
% s3 = size(x1,3);
% parfor i=1:size(x1,1)
%     for j=1:s2
%         for k=1:s3
%             D = eig([x1(i,j,k) x2(i,j,k) x3(i,j,k); x4(i,j,k) x5(i,j,k) x6(i,j,k); x7(i,j,k) x8(i,j,k) x9(i,j,k)]);
%             L1(i,j,k) = D(1);
%             L2(i,j,k) = D(2);
%             L3(i,j,k) = D(3);            
%         end
%     end
% end
%% End
end
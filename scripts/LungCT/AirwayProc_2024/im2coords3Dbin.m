function out=im2coords3Dbin(im)
% given a 3D binary matrix (1/0), return [x,y,z] coords of points with 1
im=logical(im);
n=sum(sum(sum(im)));
out=zeros(n,3); cnt=0;
for z=1:size(im,3)
    [x,y]=find(im(:,:,z));
    if ~isempty(x)
        for i=1:size(x,1)
            cnt=cnt+1;
            out(cnt,:)=[x(i) y(i) z];
        end
    end
end
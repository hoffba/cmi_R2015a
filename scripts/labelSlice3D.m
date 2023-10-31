function labelSlice3D(img,seg,voxsz,ind)

if nargin<4 || ~all(size(ind)==[1,3])
    ind = round(size(img)/2);
end

p = prctile(img(:),[0.5,99.5]);
A = cell(1,3); B = cell(1,3);
for i = 1:3
    switch i
        case 1
            a = img(:,:,ind(3));
            b = seg(:,:,ind(3));
            v = voxsz([1,2]);
        case 2
            a = squeeze(img(:,ind(2),:));
            b = squeeze(seg(:,ind(2),:));
            v = voxsz([1,3]);
        case 3
            a = squeeze(img(ind(1),:,:));
            b = squeeze(seg(ind(1),:,:));
            v = voxsz([2,3]);
    end
    a = rescale(a,0,1,InputMax=p(2),InputMin=p(1));
    A{i} = grabSlice(a,v);
    B{i} = grabSlice(labeloverlay(a,b),v);
end
figure, montage([A',B'])

function B = grabSlice(A,v)
hf = figure("Units","Pixels");
hi = imshow(A);
ha = hi.Parent;
ha.DataAspectRatio = [v,1];
ha.Units = "pixels";
pos = plotboxpos(ha);
hf.Position(3:4) = pos(3:4);
ha.Position = [0,0,pos(3:4)];
B = frame2im(getframe(ha));
delete(hf);
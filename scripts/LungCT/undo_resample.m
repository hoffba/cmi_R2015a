function undo_resample(fn_re,fn_ins)

info = niftiinfo(fn_ins);
info_re = niftiinfo(fn_re);

I = single(permute(niftiread(info_re),[2,1,3]));

d = info_re.ImageSize;
d_new = info.ImageSize;

voxsz = 0.625*ones(1,3);
voxsz_new = info.PixelDimensions;

fov = voxsz.*d;

ext = (fov - voxsz)/2;
F = griddedInterpolant({linspace(-ext(1),ext(1),d(1)),...
                        linspace(-ext(2),ext(2),d(2)),...
                        linspace(-ext(3),ext(3),d(3))},I,'nearest');

ext = (d_new-1).*voxsz_new/2;
I = F({linspace(-ext(1),ext(1),d_new(1)),...
       linspace(-ext(2),ext(2),d_new(2)),...
       linspace(-ext(3),ext(3),d_new(3))});

fn_new = erase(regexprep(fn_re,'_binVessel','.binVessel'),'.gz');
info.Datatype = 'int8';
info.BitsPerPixel = 8;
niftiwrite(int8(I),fn_new,info,'Compressed',true)
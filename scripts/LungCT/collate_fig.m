function collate_fig(hf,gifname)

[im,cm] = rgb2ind(frame2im(getframe(hf)),256);
if exist(gifname,'file')
    opts = {'WriteMode','append'};
else
    opts = {'Loopcount',Inf};
end
imwrite(im,cm,gifname,'gif',opts{:},'DelayTime',1);
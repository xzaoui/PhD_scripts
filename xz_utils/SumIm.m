
function Ims=SumIm(zim)
%%%%%%%%%%%%% sum images in a stack (or frames in a video) %%%%%%%%%%%% 
%%%% zim: stack openedwith read_stack
for i=2:length(zim)/3
    if i==2
    Ims=imadd(zim{1,1},zim{1,i});
    else
    Ims=imadd(Ims,zim{1,i},'uint16');
    end
end
    Ims=Ims-min(Ims(:));
    Ims=Ims/max(Ims(:));
imagesc(Ims);
colormap('gray')
end
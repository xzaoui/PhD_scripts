
function zim = read_stacks(imgfile,stack0,nstacks)
%%%%%%%%%%%%%%%%%%
% Name: read_stacks
% Purpose: read all the images that are in a stack or streaming video into
%           a cell object
% INPUT:
% imgfile: name of image file
% stack0: initial frame
% nstacks: total number of frames 
% OUTPUT:
% cell object with all the opened frames, cell index corrispond to frame
% index
% function developed by Alessia Lepore- El Karoui lab 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i = stack0:nstacks
       [zim{i}, map{i}] = imread(imgfile,'Index',i);
       imagesc(zim{1,1});
       [nx,ny]=size(zim{1,1});
       axis([0 ny 0 nx]);
       colormap(gray)
  end
end
%%%%%%%%%%%%%%%%%%%%  Fluorescent Single Molecule Tracking  %%%%%%%%%%%%%%%%%%%%
%%% From single-molecule movies, segment cells (by autofluorescence), 
%%% find fluorescent peaks, track them and computes Mean-Square Displacement
%%% and Diffusion Coefficient
% Designed by Xavier Zaoui, University of Edinburgh, 2015-2019
% Inspired from work by Alessia Lepore and Sebastian Jaramillo-Riveri


folder = 'movies\';         % address of movie
folderout = 'results\';     % results folder
filename='165-1509-10ms-20-1';                % name of movie


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Parameters   %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.Xtime= 150.;                                  % objective magnification                            
param.pixel2micron= 16/param.Xtime;                 % converts px in um
param.exp_time=0.010;                               % exposure time (in seconds)

%%% to segment cells
param.nstacks_track = 500;      % number of frames used to track
param.nstacks_seg = 200;        % number of frames used to segment
j0=1;                           % first frame to be read

%%% to find the peaks 
param.thrfpeak=12.5;  % multiplicative factor for trh find peak function
param.pnoise=1;       % noise in peakfind function
param.psize=4;        % size group of pixel in peakfind function
param.pgauss=7;       % size for gaussian fit peak size in centfind

%%% to track the peaks
param.maxd=8;
para.mem=0;
para.good=2;
para.dim=2;
para.quiet=0; 
param.mem=para.mem;

addpath(genpath('scriptTracking\'));        % to import additionnal functions  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Segmentation   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  for a stack (video .tiff)
imgfile = [folder filename '.tif'];
param.nstacks = param.nstacks_seg;
zim=read_stacks(imgfile,j0,param.nstacks_seg);
zs=SumIm(zim);
Im = double(mat2gray(zs));
figure, imagesc(Im), title('original image');
colormap('gray')

%%% find the threshold 
[~, threshold] = edge(Im, 'sobel');
fudgeFactor =1.0;

%%% binary gradient mask 
BWs = edge(Im,'sobel', threshold * fudgeFactor);
figure, imagesc(BWs), title('binary gradient mask');

%%% dilate the binary image
se90 = strel('disk',3); 
BWsdil = imdilate(BWs, se90);
figure, imagesc(BWsdil), title('dilated gradient mask');

%%% fill the image 
BWdfill = imfill(BWsdil, 'holes');
figure, imagesc(BWdfill)
title('binary image with filled holes');

se = strel('disk',1);
BWdfillerode = imerode(BWdfill,se);
figure, imagesc(BWdfillerode)
title('binary image with filled holes erode');

%%% clear the border 
BWnobord = imclearborder(BWdfillerode,4);
figure, imagesc(BWnobord), title('cleared border image');

%%% smoothen the obj 
seD = strel('disk',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
figure, imagesc(BWfinal), title('segmented image')

level = graythresh(Im);
bw = im2bw(BWfinal,level);
bw = bwareaopen(bw, 40 ,6);
figure, imshow(bw)

%%% outlined original image 
BWoutline = bwperim(BWfinal);
Segout = zs;
Segout(BWoutline) = 255;
figure,  imagesc(Segout),colormap(gray), title('outlined original image');

props={'perimeter','Area','Image','PixelList','Centroid','MajorAxisLength','MinorAxisLength','Orientation'};

stat_pixel=regionprops(bw,props);
stat_micron.Perimeter= stat_pixel.Perimeter*param.pixel2micron ;
stat_micron.Area=stat_pixel.Area*param.pixel2micron;
stat_micron.MajorAxisLength=stat_pixel.MajorAxisLength*param.pixel2micron;
stat_micron.MinorAxisLength=stat_pixel.MinorAxisLength*param.pixel2micron;
stat_micron.Angle=stat_pixel.Orientation;
stat_micron.Centroid=stat_pixel.Centroid;
stat_micron.PixelList=stat_pixel.PixelList;
save([folderout filename 'cellinfosMicrons.mat'],'stat_micron');

[postimeXtrack,pkpost]=findpeak_timeTOtrack(zim,param,1);
save([folderout filename 'postime2track.mat'],'postimeXtrack');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Tracking   %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totracking=postimeXtrack(:,[1:3]);
trackData_all=[];
trackData_all=track(totracking, param.maxd, para);

%%% only select tracks inside the cell area
trackData=SelTracks(trackData_all,stat_pixel);
save([folderout filename 'SELtrackData_maxdisPARA' num2str(param.maxd) '.mat'],'trackData');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Tracks Info   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

track_i=[];
for i=1:(max(trackData(:,4)))
track_i{i}= find(trackData(:,4)==i);
len_t(i)=size(track_i{i},1);  %%%%% lenght of single track , len_t>1 it's a track %%%%%%%%%%%
len_t20per(i)=size(track_i{i},1)*0.20;
end

track_info.meanlen=mean(len_t(len_t>1));
track_info.stdlen=std(len_t(len_t>1));
track_info.maxlen=max(len_t(len_t>1));
track_info.minlen=min(len_t(len_t>1));
track_info.totparticle=length(len_t);
track_info.Ntracked=length(len_t(len_t>1));
track_info.Nnotracked=length(len_t(len_t==1));
track_info.NparticleLonger4=length(len_t(len_t>=4));
track_info.meanlenLonger4=mean(len_t(len_t>=4));
track_info.stdlenLonger4=std(len_t(len_t>=4));
track_info.maxlenLonger4=max(len_t(len_t>=4));
track_info.NparticleLonger2=length(len_t(len_t>2));
track_info.meanlenLonger2=mean(len_t(len_t>2));
track_info.stdlenLonger2=std(len_t(len_t>2));
track_info.maxlenLonger2=max(len_t(len_t>2));

save([folderout filename 'track_info_'.mat'],'track_info');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Plot Tracks   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
imagesc(zim{1,1}(:,:,1))
[nx,ny]=size(zim{1,1}(:,:,1));
axis([0 ny 0 nx]);
colormap(gray)
hold on

for j=1:length(track_i)
if length(track_i{1,j})>=4
   hold all
   plot(trackData(track_i{1,j}(1:end),1),trackData(track_i{1,j}(1:end),2), 'LineWidth',1);
end
end    
savefig([folderout filename '_tracks']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot the track centered around the origin 

for j=1:length(track_i)
if length(track_i{1,j})>=4
    track2originX{i} = (trackData(track_i{1,j}(1),1)-trackData(track_i{1,j},1));
    track2originY{i} = (trackData(track_i{1,j}(1),2)-trackData(track_i{1,j},2)) ;
    track2origin{j}=cell2mat([track2originX, track2originY]);
    clear track2originX;
    clear track2originY;
end
end

figure()
for i = 1: length(track2origin)
    if length(track2origin{1,i})>1
        plot((track2origin{1,i}(:,1)),((track2origin{1,i}(:,2))),'*-','LineWidth',3);
        hold on
        xlim([-20 20]);
        ylim([-25 20]);
    end
end
savefig([folderout filename '_tracks2origin']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Pixel Distribution Occupancy   %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pix_matrix=zeros(size(zim{1,1}));

for i = 1:length(trackData)
x=round(trackData(i,1));
y=round(trackData(i,2)); 
if  x>0 && x <=size(zim{1,1},2) && y>0 && y <=size(zim{1,1},1)
    if pix_matrix(y,x)==0
       pix_matrix(y,x)= +1 ;
    else
       pix_matrix(y,x)= pix_matrix(y,x)+1;

    end
end
end

figure();
title('position distribution within pixels');
imagesc(pix_matrix);
[nx,ny]=size(zim{1,1});
axis([0 ny 0 nx]);
colorbar;
savefig([folderout filename 'pos_2Ddist']); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Square Displacement & Displacement   %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sd_all,sd_micron]=SDcalculations(trackData,param); 
save([folderout filename 'SDmicron.mat'],'sd_micron'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mean Square Displacement %%%

[msd,msd_micron]=MSDcalculations(sd_all,sd_micron); 
save([folderout filename 'MSDmicron.mat'],'msd_micron');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot MSD n time lag +errors bars 

msdM=msd;
msd_micronM1=msd_micron;

for i =1:1:size(msd_micronM1,1)
       msd_micronM2(i,:)=msd_micronM1(i,1:(end-1));       
end

figure()
for i =1:1:size(msd_micronM1,1)
   hold all
   yi=msd_micronM2(i,1:end);
   timex=(1:length(yi(yi>0)))*param.exp_time*1000;
    plot( timex,yi(yi>0),'*-','LineWidth',2);
    xlabel('time(ms)')
    ylabel('MSD (micron^2)')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Diffusion Coefficient  %%%

D_1step=D_1t_step_calculation(msd_micron,param);
save([folderout filename 'D_1timestep.mat'],'D_1step');    

figure()
histogram(log10(D_1step),'BinWidth', 0.10, 'Normalization','probability')
ylabel('frequency','FontSize',20)
xlabel('log(D_{app}) (\mum^{2}s^{-1})','FontSize',20)
title(['N=' num2str(length(D_1step)) ';' filename ])
savefig([folderout filename 'maxdisp_PARA' num2str(param.maxd) 'D_1timestepLOGHistogram']);


Dbins = [logspace(-2,2,20)];  % edges of the histogram count: logspace(A,B,N) =  N increments between 10^A and 10^B
av_D = mean(D_1step);

figure()
h = plot(Dbins,relhist(D_1step,Dbins),'LineWidth',3);
axis([1e-2 1e2 0 0.5]);
set(gca,'xscale','log');
set(gca,'XTick',[1e-1,1e0,1e01,1e02]);
set(gca,'YTick',[0,0.1, 0.2, 0.3, 0.4, 0.5]);
ylabel('frequency','FontSize',20)
xlabel('D_{app} (\mum^{2}s^{-1})','FontSize',20)    
title(['N=' num2str(length(D_1step)) ';' filename ])
a = annotation('textbox',[.15 .8 .1 .1],'String',['mean = ' num2str(av_D)],'FitBoxToText','on','Color','blue');
savefig([folderout filename 'maxdisp_PARA' num2str(param.maxd) 'D_1timesLOGplot']);

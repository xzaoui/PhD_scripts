 function [postimeXtrack,pos]=findpeak_timeTotrack(zim,param,doplot)
%%%%%%%%%%%%%%%%%%
% Name: findpeak_timeTotrack
% Purpose: find the peaks positions and prepare the file to give to the
% track function
%           
% INPUT:
% zim: images organised in a cell file (see read_stacks)
% param: parameters for find the peaks 
%    param.thrfpeak=12.5;  %% multiplicative factor for trh find peak function
%    param.pnoise=1;       %% noise in peakfind function
%    param.psize=4;        %% size group of pixel in peakfind function
%    param.pgauss=7;       %% size for gaussian fit peak size in centfind 
% OUTPUT:
% file(table) organize like
%                    (x)      (y)      (t)
% ;     pos = 3.60000      5.00000      0.00000
% ;           15.1000      22.6000      0.00000
% ;           4.10000      5.50000      1.00000 
% ;           15.9000      20.7000      2.00000
% ;           6.20000      4.30000      2.00000
% to give to the track function
%
% function developed by Alessia Lepore- El Karoui lab 2017
% other functions used here(see documentation): 
% pkfind
% bpass
% centfind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% filter image and find peaks %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imfilt=[];
    pkpos=[];
    for i = 1 : param.nstacks
       % thr=mean(mean(zim{i,1}));
        %thr=std(mean(imfilt{i,1}));
        imfilt{i}=bpass(zim{1,i},param.pnoise,param.psize);             
        thr=mean(mean(imfilt{i}));
        pkpos{i}=pkfnd(imfilt{i},thr*param.thrfpeak,param.psize);    
   
        if length(pkpos{i}) > 1 
           [pkcentCoords{i},pkcentCoords_extOut{i}]=centfind(imfilt{i},(pkpos{i}),param.pgauss,0,['Gaussian','interactive']);

        end
    end

    pkpost=(pkpos);
    for i =1 : length (pkcentCoords)
        if length(pkcentCoords{i})>1
            pkpost{i}= double(pkcentCoords{i}(:,1:2));
        end
    end
   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% make plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if doplot==1
        
    figure()
    for i = 1 :50:(param.nstacks)
           if length(pkpost{1,i})>1
            if i==1
                subplot(ceil(param.nstacks/(50*2)),2,i)
%                 imagesc(zim{1,i})
%                 colormap('gray')
%                [ny,nx]=size(zim{1,1});
                imagesc( imfilt{i})
                colormap('gray')
                 [ny,nx]=size( imfilt{i});
    %             imagesc( zim_crop1{i,1})
    %             colormap('gray')
    %             [ny,nx]=size( zim_crop1{i,1});
    %             axis([0 nx 0 ny]);
                hold all
                plot(pkpost{1,i}(:,1),pkpost{1,i}(:,2),'.r','MarkerSize',8)
    %             plot(pkcentCoords{1,i}(:,1),pkcentCoords{1,i}(:,2),'.','MarkerSize',8)
                figure()
                subplot(2,ceil(param.nstacks/(50*2)),i)
            else
                subplot(ceil(param.nstacks/(50*2)),2,fix(i/50))
%                 imagesc(zim{1,i})
%                 colormap('gray')
%                [ny,nx]=size(zim{1,1});
                 imagesc( imfilt{i})
                 colormap('gray')
                 [ny,nx]=size( imfilt{i});
    %             imagesc( zim_crop1{i})
    %             colormap('gray')
    %             [ny,nx]=size( zim_crop1{i});
    %             axis([0 nx 0 ny]);
                hold all
                plot(pkpost{1,i}(:,1),pkpost{1,i}(:,2),'.r','MarkerSize',8)
    %             plot( pkcentCoords{1,i}(:,1),pkcentCoords{1,i}(:,2),'.','MarkerSize',8)
            end    
        end
    end
    end
    
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% make files with peak pos and time to track%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
    
    pos=transpose(pkpost);
    posxtrack=cell2mat(transpose(pkpost));
    % pos=(pkcentCoords);
    % posxtrack=cell2mat((pkcentCoords));


    for i=1:param.nstacks
        if length(pkpos{1,i})>1
        time{i}=[ones(size(pkpost{1,i},1),2)*i*param.exp_time];
        end
    end

    ttime=transpose(time);
    postime= [pos ttime];
    postimeXtrack = cell2mat(postime);
    postimeXtrack(:,4)=0;

     
    
    
    end
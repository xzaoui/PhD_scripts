function [sd_all,sd_micron]=SDcalculations(trackData,param)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Name: SDcalculation in 2D
    %Purpoise:compute the square displacement from tracked data
    %method:  square displacement (sd) for all the n time intervals only
    %         for object tracked for more than 4 frames
    %INPUT: 
    %trackData (usually positions and time processed with function track)
    %trackData format
    %    (x)            (y)         (t)        (id_object)
    %    3.60000      5.00000      0.00000      0.00000
    %    4.10000      5.50000      1.00000      0.00000
    %    6.20000      4.30000      2.00000      0.00000
    %    15.1000      22.6000      0.00000      1.00000
    %    15.9000      20.7000      2.00000      1.00000
    %param :  param.pixel2micron conversion factor from pixel to micron
    %         param.pixel_size/param.Xtime(objective magnification);
    %OUTPUT:
    %sd_all : square displacement for all the possible time intervals of
    %         trajetories longer than 4 frames
    %sd_micron: same of sd_all in micron
    %function developed by Alessia Lepore- El Karoui lab 2017
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    track_i=[];
    for i=1:(max(trackData(:,4)))
        track_i{i}= find(trackData(:,4)==i);
        len_t(i)=size(track_i{i},1);  %%%%% lenght of single track , len_t>1 it's a track %%%%%%%%%%%
    end
    
    k=1;
    track_ind=[];
    for i=1:length(track_i)
       if (length(track_i{1,i}))>=4.
           track_ind{1,k}=track_i{1,i};  
           k=1+k;
       end
    end
    
    for i=1:length(track_ind)
       if (length(track_ind{1,i}))>=4.
          for n = 1:length(track_ind{1,i})
           sd=[];
           for j = 1: (length(track_ind{1,i})-n)    
               sd(j)=((trackData(track_ind{1,i}(j+n),1)-trackData(track_ind{1,i}(j),1)).^2+(trackData(track_ind{1,i}(j+n),2)-trackData(track_ind{1,i}(j),2)).^2); 
           end
            sd_all{i,n}=sd;
            sd_micron{i,n}=sd*param.pixel2micron.^2;
          end
        
        end
    end
  
    
    
end
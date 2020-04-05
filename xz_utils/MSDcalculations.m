function [msd,msd_micron]=MSDcalculations(sd_all,sd_micron)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Name: MSDcalculation in 2D
    %Purpoise:compute the mean square displacement from tracked data
    %method:  mean square displacement (msd) for all the n time intervals only
    %         for object tracked for more than 4 frames using the sd
    %         computed by SDcalculations function
    %INPUT: 
    %sd_all : square displacement for all the possible time intervals of
    %         trajetories longer than 4 frames
    %sd_micron: same of sd_all in micron
    %
    %OUTPUT:
    % msd: mean square displacement for track longer than 4 step 
    %msd_micron: mean square displacement for track longer than 4 step in
    %            micron
    %
    %function developed by Alessia Lepore- El Karoui lab 2017
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    msd=[];
    %msd_std=[];
    msd_micron=[];
    %msd_micronstd=[];
    
    for i=1:size(sd_all,1)
       for n = 1: size(sd_all,2)
           if sd_all{i,n}>0
             msd(i,n)=mean(sd_all{i,n});
             %msd_std(i,n)=std(sd_all{i,n});
             msd_micron(i,n)=mean(sd_micron{i,n});
             %msd_micronstd(i,n)=std(sd_micron{i,n});
           end
       end

    end
end
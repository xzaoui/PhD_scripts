function trackIN=SelTracks(trackData,stat_pixel)
%%%%% select the track only inside the cell area %%%%%%%%%%
%%% INPUT 
%%% trackData: all the tracked spot in the format
%%%       (x)            (y)         (t)        (id_object)
%    %    3.60000      5.00000      0.00000      0.00000
%    %    4.10000      5.50000      1.00000      0.00000
%    %    6.20000      4.30000      2.00000      0.00000
%    %    15.1000      22.6000      0.00000      1.00000
%    %    15.9000      20.7000      2.00000      1.00000
%%% stat_pixel contain info about the segmented cell, in particular PixelList
%%%(before use this faction segment!)
%%% OUTPUT
%%% trackIN: track of objects inside the cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 j=1
  for i=1:size(trackData,1) 
       if ismember(ceil(trackData(i,1)),stat_pixel.PixelList(:,1)) == 1 && ismember(ceil(trackData(i,2)),stat_pixel.PixelList(:,2))==1;
           trackIN(j, 1:4)= trackData(i,1:4);
           j=j+1;
       end
  end

end
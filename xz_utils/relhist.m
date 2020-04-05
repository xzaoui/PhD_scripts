function [hf] = relhist(x,bins)   % normalise counts by total number

h = histc(x,bins);
hf = h./sum(h);

end
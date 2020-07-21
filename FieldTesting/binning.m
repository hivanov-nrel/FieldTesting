function [bin_center,var_count,var_average,var_std] = binning(bin_chan, bin_bounds, var)
%Jason Roadman
%NREL NWTC
%6/10/13
%This function conducts simple binning analysis of a variable (var) against the bin channel (bin_chan) over
%the bin spacing (bin_bounds)
%INPUTS:    bin_chan - vector of channel to bin against (typically wind speed)
%           bin_bounds - vector of boundaries for bins
%               note: two additional bins will be created, less than or equal to min(bin_bounds) and greater than max(bin_bounds)
%           var - vector of data to be averaged in each bin
%               (if not explicitly specified, will be the same as bin_chan, i.e. bin WS against itself.)
%OUTPUTS:   var_count - number of elements of var in each bin
%           var_average - bin averages of var
%USAGE:
% [var_count,var_average, bin_center] = binning(bin_chan, bin_bounds, var)
% [var_count,var_average, bin_center] = binning(bin_chan, bin_bounds)


%bin bin_chan if no other variable is specified
if nargin <3
    var = bin_chan;
end

%find elements less than or equal to first bin
i_less = find(bin_chan <= bin_bounds(1));
if isempty(i_less)
    var_count(1,1) = 0;
    var_average(1,1) = 0;
    var_std(1,1) = 0;
else
    var_count(1,1) = length(i_less);
    var_average(1,1) = mean(var(i_less));
    var_std(1,1) = std(var(i_less),1);
end
bin_center(1,1) = -inf;

%loop through each bin
for ii = 1:length(bin_bounds)-1
   i_bin = find(bin_bounds(ii) < bin_chan & bin_chan <= bin_bounds(ii+1));
   bin_center(ii+1,1) = (bin_bounds(ii)+bin_bounds(ii+1))/2;
   if isempty(i_bin)
       var_count(ii+1,1) = 0;
       var_average(ii+1,1) = 0;
       var_std(ii+1,1) = 0;
   else
       var_count(ii+1,1) = length(i_bin);
       var_average(ii+1,1) = mean(var(i_bin));
       var_std(ii+1,1) = std(var(i_bin),1);
   end
end

%find elements greater than last bin
i_more = find(bin_chan > bin_bounds(end));
if isempty(i_more)
    var_count(end+1,1) = 0;
    var_average(end+1,1) = 0;
    var_std(end+1,1) = 0;
else
    var_count(end+1,1) = length(i_more);
    var_average(end+1,1) = mean(var(i_more));
    var_std(end+1,1) = std(var(i_more),1);
end
bin_center(end+1,1) = inf;

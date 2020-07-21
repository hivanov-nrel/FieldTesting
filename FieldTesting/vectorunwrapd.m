function [data_unwrapped] = vectorunwrapd(data)
%vectorunwrap
%Jason Roadman
%NREL NWTC
%11/17/11
%This function takes a vector of numbers, presumably a directional signal in degrees
% and ensures that every value is between 0 and 360.

for ii = 1:length(data)
	if(data(ii)<0)
		data(ii) = data(ii)+360;
	elseif(data(ii)>360)
		data(ii) = data(ii)-360;
    end
end

if(max(data) > 360 | min(data) < 0)
    data = vectorunwrapd(data);
end

data_unwrapped = data;

function [captmat,TIbinbound,WSbinbound] = gencaptmat(stats,windspeed,i_desired)
%INPUTS:
%   stats - statistics data structure in the standard field test format
%   windspeed - string corresponding to the wind speed channel of interest
%                   i.e. "Controller_Nacelle_Ano"
%   i_desired - indices of interest (i.e. so one can filter stats for online)
%               *optional input
%OUTPUTS:
%   captmat - capture matrix containing number of stats in each bin
%               NOTE: this matrix will not have headers
%   TIbinbound - boundaries of turbulence intensity (i.e. rows of captmat)
%   WSbinbound - boundaries of wind speed (i.e. header of captmat)
%TODO:
%   - add functionality to input different bin boundaries as optional inputs
%   - add functionality to save off directly as an excel or other type of data 
%           *if not within this function, then a better companion function than the patch written during the Alstom foundation test



%constants
WSbinbound = 2:1:25; %(2,3],(3,4],...,(24:25],(25,inf)
TIbinbound = [3:2:29]'; %(-inf,3],(3,5],(5,7],...,(27:29],(29,inf)



%initialize
captmat = zeros(length(TIbinbound)+1,length(WSbinbound));

if nargin < 3
    i_desired = 1:length(stats.(windspeed).mean10min);
end

%calculate turbulence intensity
TI = stats.(windspeed).stdev10min./stats.(windspeed).mean10min*100;



%% generate main portion of capture matrix (minus first row, last row, and last column
%loop through each WS bin
for fw = 1:length(WSbinbound)-1
    
    %loop through each TI bin
    for ft = 1:length(TIbinbound)-1
                
        %find subset that lives in WS bin
        i_WS = find(WSbinbound(fw) < stats.(windspeed).mean10min(i_desired) & stats.(windspeed).mean10min(i_desired) <= WSbinbound(fw+1));

        %find subset that lives in TI bin
        i_TI = find(TIbinbound(ft) < TI(i_desired) & TI(i_desired) <= TIbinbound(ft+1));
        
        %calculate number that exist in overall bin and populate capture matrix accordingly
        captmat(ft+1,fw) = length(intersect(i_WS,i_TI));        

    end

end

%% generate first row of capture matrix minus last cell

%loop through each WS bin
for fw = 1:length(WSbinbound)-1

        %find subset that lives in WS bin
        i_WS = find(WSbinbound(fw) < stats.(windspeed).mean10min(i_desired) & stats.(windspeed).mean10min(i_desired) <= WSbinbound(fw+1));

        %find subset that lives in TI bin
        i_TI = find(TI(i_desired) <= TIbinbound(1));
        
        %calculate number that exist in overall bin and populate capture matrix accordingly
        captmat(1,fw) = length(intersect(i_WS,i_TI));        

end

%% generate last row of capture matrix minus last cell
for fw = 1:length(WSbinbound)-1

        %find subset that lives in WS bin
        i_WS = find(WSbinbound(fw) < stats.(windspeed).mean10min(i_desired) & stats.(windspeed).mean10min(i_desired) <= WSbinbound(fw+1));

        %find subset that lives in TI bin
        i_TI = find(TIbinbound(end) < TI(i_desired));
        
        %calculate number that exist in overall bin and populate capture matrix accordingly
        captmat(length(TIbinbound)+1,fw) = length(intersect(i_WS,i_TI));        

end

%% generate last column of capture matrix minus first and last cell

%loop through each TI bin
for ft = 1:length(TIbinbound)-1

    %find subset that lives in WS bin
    i_WS = find(WSbinbound(end) < stats.(windspeed).mean10min(i_desired));

    %find subset that lives in TI bin
    i_TI = find(TIbinbound(ft) < TI(i_desired) & TI(i_desired) <= TIbinbound(ft+1));

    %calculate number that exist in overall bin and populate capture matrix accordingly
    captmat(ft+1,end) = length(intersect(i_WS,i_TI));        
end


%% generate first and last cells of last column of capture matrix

%upper righthand cell
i_WS = find(WSbinbound(end) < stats.(windspeed).mean10min(i_desired));
i_TI = find(TI(i_desired) <= TIbinbound(1));
captmat(1,end) = length(intersect(i_WS,i_TI));    

%lower righthand cell
i_WS = find(WSbinbound(end) < stats.(windspeed).mean10min(i_desired));
i_TI = find(TIbinbound(end) < TI(i_desired));
captmat(end,end) = length(intersect(i_WS,i_TI));    



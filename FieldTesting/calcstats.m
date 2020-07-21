function [data] = calcstats(data,statlength,statslist)
%Jason Roadman
%NREL NWTC
%9/4/12
%This function takes a Matlab data file (in the standard format for the Field Testing Group)
%  and calculates (or recalculates and overwrites) statistics for each processed channel.  
%  Optionally, the user can specify the list of channels they would like statistics calculated on.
%INPUT:     data - struct in the standard format for the Field Testing Group
%           statlength - string specifying averaging time, 
%                           *must be one of three values: '10min', '1min', '10sec'
%           statlist - optional third input, class array of channels that need statistics calculated on
%OUTPUT:    data - same struct, with statistics fields updated
%USAGE:     [data] = calcstats(data,statlength): recalculate all channels
%           [data] = calcstats(data,statlength,statslist): recalculate stats only for those listed in statslist
%WARNINGS:  
%     - ***DANGER: if a sample is dropped (scan error, etc) the averages will incorporate this skip i.e. 1:00, 2:01, ... instead of 1:00, 2:00, ....
%     - the code assumes all channels are sampled at the same rate and are the same length
%     - the code assumes all channels requiring vector averages are in degrees.  It returns their angle relative to 
%           the compass rose (i.e. 0 deg = north, 90 deg = east, etc)
%     - the code also assumes that all vector channels have been unwrapped into the 0-360 deg range.  Results could be wrong if not.
%             - This was achieved as a defualt step in INITFTMAT after version 179
%TODO:
%     - fix the error described above in the first warning by actually finding averaging periods based on the time stamp, not number of samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%default case: process all channels
if nargin < 3
    statslist = data.chanlist;
end

%flag for determining if file is long enough to calculate desired statistical period on
lengthflag = false;

logentry = ['AC - Stats calculated using "' mfilename '.m" - ' datestr(now)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%determine which type of averaging period to implement
if isequal(statlength, '10min')
    %determine averaging period, in seconds
    avperiod = 10*60; %seconds, 10 min = 600 sec
%     warning('10 minute statistics not tested with new version yet.')
elseif isequal(statlength, '1min')
    %determine averaging period, in seconds
    avperiod = 1*60; %seconds, 1 min = 60 sec
%     warning('1 minute statistics not tested yet.')
elseif isequal(statlength, '10sec')
    %determine averaging period, in seconds
    avperiod = 10; %seconds
%     warning('10 second statistics not tested yet.')
else
    error('%s is not a valid input for averaging period.')
end



%name fields based on averaging type requested
meanfield = ['mean' statlength];
maxfield = ['max' statlength];
minfield = ['min' statlength];
stdevfield = ['stdev' statlength];


%calculate number of points in averaging sample
avsample = double(avperiod*50);

%check if file is long enough to perform average, 
%set flag to initialize to NaN values if it's not
%calculate indices of interest for averages if it is
if length(data.(data.chanlist{end}).scaled_data) < avsample %%% EDIT: using scaled_data
    warning('File is not long enough to calculate average desired, NaN stored in stats.')
    
    %set flag to initialize to NaN for stats
    lengthflag = true;
    
else 
    %calculate length of processed array
    numpts = length(data.(data.chanlist{10}).scaled_data); %%% EDIT: using scaled_data
    
    %Build array of indices where each row contains all the indices for a single average, 
    %   total number of rows equal to the number of averages desired
    i_av = [1:avsample:floor(numpts/avsample)*avsample]'*ones(1,avsample)+ones(floor(numpts/avsample),1)*[0:avsample-1];
    %fprintf('Calculating %s statistics.\n',statlength)
end


%loop through each channel in statslist for 10 minute stats 
%%%WARNING: this code only has been checked for 10 minute stats
for ii = 1:length(statslist)

    %if file is too short to calculate desired average period, store NaN in stats
    if lengthflag 
          
        %average
        data.(statslist{ii}) = setfield(data.(statslist{ii}),meanfield, NaN);
            
        %max
        data.(statslist{ii}) = setfield(data.(statslist{ii}),maxfield, NaN);
            
        %min
        data.(statslist{ii}) = setfield(data.(statslist{ii}),minfield, NaN);
        
        %standard deviation
        data.(statslist{ii}) = setfield(data.(statslist{ii}),stdevfield, NaN);
        
        %add to log field
        %data.(statslist{ii}).log{end+1,1} = ['AD - File is not long enough to calculate average desired, NaN stored in stats, ' datestr(now)];
        
        %skip to next loop iteration
        continue
    end
        
    %grab data chunk --> easier to code
    datachunk = data.(statslist{ii}).scaled_data(i_av); %%% EDIT: read scaled_data instead of data
    
    %maks sure datachunk is a row vector
    if size(datachunk,1) ~= 1 && size(datachunk,2) == 1
        datachunk = datachunk';
    end
        
    %check if channel is a time channel
    if isequal(data.(statslist{ii}).avtype, 'time')
%         fprintf('%s == time\n',statslist{ii})
        
        %average is the first value from the time period
        data.(statslist{ii}) = setfield(data.(statslist{ii}),meanfield, datachunk(:,1));
            
        %max
        data.(statslist{ii}) = setfield(data.(statslist{ii}),maxfield, max(datachunk,[],2));
            
        %min
        data.(statslist{ii}) = setfield(data.(statslist{ii}),minfield, min(datachunk,[],2));
        
        %standard deviation
        data.(statslist{ii}) = setfield(data.(statslist{ii}),stdevfield, std(datachunk,0,2));
        
        %add to log field
        %data.(statslist{ii}).log{end+1,1} = logentry;

    %check if channel is a vector averaging channel
    elseif isequal(data.(statslist{ii}).avtype, 'vector')
%             fprintf('%s == vector\n',statslist{ii})
        




        %allocate memory for stdev and mean
        data.(statslist{ii}) = setfield(data.(statslist{ii}),meanfield,-ones(size(datachunk,1),1));
        data.(statslist{ii}) = setfield(data.(statslist{ii}),stdevfield,-ones(size(datachunk,1),1));
        
        %loop through each averaging period to calculate mean and stdev 
        for aa = 1:size(datachunk,1)
            
            %grab the appropriate row of datachunk to calculate a single average and standard deviation on
            smallchunk = datachunk(aa,:);
        

            %mean (use routine from Cambell data logger)
            Ux = sum(sind(smallchunk))/length(smallchunk);
            Uy = sum(cosd(smallchunk))/length(smallchunk);
            vector_av = 90 - atan2(Uy,Ux)*180/pi();
            if(vector_av < 0)
                vector_av = vector_av+360;
            elseif(vector_av>360)
                vector_av = vector_av-360;
            end
            data.(statslist{ii}).(meanfield)(aa) = vector_av;
        
            
            %standard deviation (use routine from Cambell data logger/Yamartino algorithm)
            magsum = round((Ux^2 + Uy^2)*1e8)/1e8; %round to the 8th decimal place (10 ppb) to help reduce roundoff error
            epsilon = ( 1 - magsum)^0.5;

            if  ~isreal(epsilon) 
                vector_stdev = 0;
                warnmsg = ['epsilon less than zero in Yamartino algorithm. "' statslist{ii} '" assumed constant and stdev set to zero.']

                %add to log field
                data.(statslist{ii}).log{end+1,1} = ['AE - (Ux^2+Uy^2) = ', num2str(magsum,'%12.10f'), ' in Yamartino algorithm. Signal assumed constant and stdev set to zero. ', datestr(now)];

                warning(warnmsg)

                fprintf('(Ux^2+Uy^2) = %12.10f\n',(Ux^2+Uy^2))

            else
                vector_stdev = asind(epsilon)*(1+0.1547*epsilon^3);
            end
            data.(statslist{ii}).(stdevfield)(aa) = vector_stdev;

        end
        
        
        
        
        
        
                
        %max
        data.(statslist{ii}) = setfield(data.(statslist{ii}),maxfield, max(datachunk,[],2));
            
        %min
        data.(statslist{ii}) = setfield(data.(statslist{ii}),minfield, min(datachunk,[],2));
        
        %add to log field
        %data.(statslist{ii}).log{end+1,1} = logentry;
        
    %check if channel is a scalar averaging channel
    elseif isequal(data.(statslist{ii}).avtype, 'scalar')
%             fprintf('%s == scalar\n',statslist{ii})

        %average
        data.(statslist{ii}) = setfield(data.(statslist{ii}),meanfield, mean(datachunk,2));
            
        %max
        data.(statslist{ii}) = setfield(data.(statslist{ii}),maxfield, max(datachunk,[],2));
            
        %min
        data.(statslist{ii}) = setfield(data.(statslist{ii}),minfield, min(datachunk,[],2));
        
        %standard deviation
        data.(statslist{ii}) = setfield(data.(statslist{ii}),stdevfield, std(datachunk,0,2));
        
        %add to log field
        %data.(statslist{ii}).log{end+1,1} = logentry;


    %error if one of three above options is not met
    else
        errmsg = ['Error: channel "' statslist{ii} '" not classified correctly for averaging.'];
        error(errmsg)
    end

    
    clear datachunk
    
end

    
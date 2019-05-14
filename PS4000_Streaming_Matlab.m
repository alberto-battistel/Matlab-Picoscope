%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Filename:    PS4000_IC_Generic_Driver_Streaming
%
% Copyright:   Pico Technology Limited 2014
%
% Author:      KPV
%
% Description:
%   This is a MATLAB script that demonstrates how to use the
%   PicoScope 4000 Series Instrument Control Toobox driver to 
%   collect data in streaming mode for 2 channels without aggregation or a 
%   trigger.
%
%	To run this application:
%		Ensure that the following files/folders are located either in the 
%       same directory or define the path in the PS4000Config.m file:
%       
%       - picotech_ps4000_generic.mdd
%       - PS4000Constants
%       - ps4000.dll & ps4000Wrap.dll 
%       - ps4000MFile.m & ps4000WrapMFile.m
%       - PicoConstants.m
%       - PicoStatus.m
%       - Functions
%
%   Device used to generated example: PicoScope 4423 & 4262
%
% modified several times by Alberto Battistel
% last time on 01 April 2016
% add all the folder pico_V2 with its own subfolders to matlab path
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Workspace

clear;
clc;
% cwd = pwd;
% cd(cwd)
% addpath('..')

%% Alberto
Tinterval = 2.5e-6;       % s
Nsample = 10e6;      % number points
folder2save = pwd;
% change voltage scale on line 80 something

%% Load Configuration Information

PS4000Config;

%% Parameter Definitions
% Define any parameters that might be required throughout the script.

channelA = ps4000Enuminfo.enPS4000Channel.PS4000_CHANNEL_A;
channelB = ps4000Enuminfo.enPS4000Channel.PS4000_CHANNEL_B;

%% Device Connection

% Create device -  specify serial number if required
% Specify serial number as 2nd argument if required.
ps4000DeviceObj = icdevice('picotech_ps4000_generic' ); 
% 'AW737/010'
% Connect device
connect(ps4000DeviceObj);

%% Display Unit Information

[infoStatus, unitInfo] = invoke(ps4000DeviceObj, 'getUnitInfo');





%% Channel Setup
% All channels are enabled by default - switch off all except Channels A and B.

numChannels = get(ps4000DeviceObj, 'channelCount');

% Channel A
channelSettings(1).enabled = PicoConstants.TRUE;
channelSettings(1).coupling = 1;
channelSettings(1).range = ps4000Enuminfo.enPS4000Range.PS4000_2V;      % Alberto A range

channelARangeMV = PicoConstants.SCOPE_INPUT_RANGES(channelSettings(1).range + 1);

% Channel B
channelSettings(2).enabled = PicoConstants.TRUE;
channelSettings(2).coupling = 1;
channelSettings(2).range = ps4000Enuminfo.enPS4000Range.PS4000_500MV;     % Alberto B range

channelBRangeMV = PicoConstants.SCOPE_INPUT_RANGES(channelSettings(2).range + 1);


if (numChannels == 4)

    % Channel C
    channelSettings(3).enabled = PicoConstants.FALSE;
    channelSettings(3).coupling = 1;
    channelSettings(3).range = ps4000Enuminfo.enPS4000Range.PS4000_2V;


    % Channel D
    channelSettings(4).enabled = PicoConstants.FALSE;
    channelSettings(4).coupling = 1;
    channelSettings(4).range = ps4000Enuminfo.enPS4000Range.PS4000_2V;
    

end

% Keep the status values returned from the driver.
for ch = 1:numChannels
   
    status.setChannelStatus(ch) = invoke(ps4000DeviceObj, 'ps4000SetChannel', ...
        (ch - 1), channelSettings(ch).enabled, ...
        channelSettings(ch).coupling, channelSettings(ch).range);
    
end

% Obtain the maximum ADC Count from the driver
maxADCCount = double(get(ps4000DeviceObj, 'maxADCValue'));

%% Trigger Setup
% Turn off trigger
% If a trigger is set and the autoStop property in the driver is set to
% '1', the device will stop collecting data once the number of post trigger
% samples have been collected.

% Trigger properties and functions are located in the Instrument
% Driver's Trigger group.

triggerGroupObj = get(ps4000DeviceObj, 'Trigger');
triggerGroupObj = triggerGroupObj(1);

[status.setTriggerOff] = invoke(triggerGroupObj, 'setTriggerOff');

%% Set Data Buffers
% Data buffers for Channel A and B - buffers should be set with the driver,
% and these MUST be passed with application buffers to the wrapper driver
% in order to ensure data is correctly copied.

overviewBufferSize  = 250000; % Size of the buffer to collect data from buffer.
segmentIndex        = 0;   

ratioMode = ps4000Enuminfo.enRatioMode.RATIO_MODE_NONE;

% Buffers to be passed to the driver
pDriverBufferChA = libpointer('int16Ptr', zeros(overviewBufferSize, 1));
pDriverBufferChB = libpointer('int16Ptr', zeros(overviewBufferSize, 1));

status.setDataBufferChA = invoke(ps4000DeviceObj, 'ps4000SetDataBufferWithMode', ...
    channelA, pDriverBufferChA, overviewBufferSize, segmentIndex, ratioMode);

status.setDataBufferChB = invoke(ps4000DeviceObj, 'ps4000SetDataBufferWithMode', ...
   channelB, pDriverBufferChB, overviewBufferSize, segmentIndex, ratioMode);

% Application Buffers - these are for copying from the driver into temporarily.
pAppBufferChA = libpointer('int16Ptr', zeros(overviewBufferSize, 1));
pAppBufferChB = libpointer('int16Ptr', zeros(overviewBufferSize, 1));

% Streaming properties and functions are located in the Instrument
% Driver's Streaming group.

streamingGroupObj = get(ps4000DeviceObj, 'Streaming');
streamingGroupObj = streamingGroupObj(1);

% Register application buffer and driver buffers (with the wrapper).

status.setAppAndDriverBuffersA = invoke(streamingGroupObj, 'setAppAndDriverBuffers', channelA, ...
    pAppBufferChA, pDriverBufferChA, overviewBufferSize);

status.setAppAndDriverBuffersB = invoke(streamingGroupObj, 'setAppAndDriverBuffers', channelB, ...
   pAppBufferChB, pDriverBufferChB, overviewBufferSize);

%% Run Streaming And Get Values
% Use default value for streaming interval which is 1e-6 for 1MS/s
% Collect data for 1 second with auto stop - maximum array size will depend
% on PC's resources - type 'memory' at MATLAB command prompt for further
% information.

% To change the sample interval e.g 5 us for 200KS/s. The call to
% ps4000RunStreaming will output the actual sampling interval used by the 
% driver.
% Tinterval = 2e-6;
set(streamingGroupObj, 'streamingInterval', Tinterval);

% For 10MS/s, specify 100ns
%set(streamingGroupObj, 'streamingInterval', 100e-9);

% Set the number of pre- and post-trigger samples
% If no trigger is set 'numPreTriggerSamples' is ignored
% Nsample = 80e6;      % Alberto
set(ps4000DeviceObj, 'numPreTriggerSamples', 0);
set(ps4000DeviceObj, 'numPostTriggerSamples', Nsample);     % Alberto

% Alberto
fprintf('\n\ntime: %.1d min\n\n\n', Tinterval*Nsample/60)
fprintf('\n\nNsamples: %.3g\n\n\n', Nsample)
fprintf('\n\nTinterval: %.3g s\n\n\n', Tinterval)

% autoStop parameter can be set to false (0)
%set(streamingGroupObj, 'autoStop', PicoConstants.FALSE);

% Set other streaming parameters
downSampleRatio     = 1;

% Defined buffers to store data collected from channels.
% If capturing data without using the autoStop flag, or if using a trigger 
% with the autoStop flag, allocate sufficient space (1.5 times the size is 
% shown below) to allow for pre-trigger data. Pre-allocating the array is 
% more efficient than using vertcat to combine data.

maxSamples = get(ps4000DeviceObj, 'numPreTriggerSamples') + ...
    get(ps4000DeviceObj, 'numPostTriggerSamples');

% Take into account the downSamplesRatioMode - required if collecting data
% without a trigger and using the autoStop flag.
finalBufferLength = round(1.5 * maxSamples / downSampleRatio);

pBufferChAFinal = libpointer('int16Ptr', zeros(maxSamples, 1, 'int16'));
pBufferChBFinal = libpointer('int16Ptr', zeros(maxSamples, 1, 'int16'));

% % Alberto
% A = zeros(maxSamples, 1, 'int16');
% B = zeros(maxSamples, 1, 'int16');
nBuffer = ceil(maxSamples/overviewBufferSize);  % try to preallocare a cell
bufferChAmV = cell(10e3,1);
bufferChBmV = cell(10e3,1);

% Prompt to press enter to begin capture
input('Press ENTER to begin data collection.', 's');

[status.runStreaming, actualSampleInterval, sampleIntervalTimeUnitsStr] = ...
    invoke(streamingGroupObj, 'ps4000RunStreaming', downSampleRatio, overviewBufferSize);
    
disp('Streaming data...');
fprintf('Click the STOP button to stop capture or wait for auto stop if enabled.\n\n') 

% Variables to be used when collecting the data
hasAutoStopped      = PicoConstants.FALSE;
newSamples          = 0; % Number of new samples returned from the driver.
previousTotal       = 0; % The previous total number of samples.
totalSamples        = 0; % Total samples captured by the device.
startIndex          = 0; % Start index of data in the buffer returned.
hasTriggered        = 0; % To indicate if trigger has occurred.
triggeredAtIndex    = 0; % The index in the overall buffer where the trigger occurred.

time = zeros(overviewBufferSize / downSampleRatio, 1);	% Array to hold time values

status.getStreamingLatestValues = PicoStatus.PICO_OK; % OK

[stopFig.h, stopFig.h] = stopButton();             
             
flag = 1; % Use flag variable to indicate if stop button has been clicked (0)
setappdata(gcf, 'run', flag);

% % Plot Properties
% 
% % Plot on a single figure
% figure1 = figure;
% axes1 = axes('Parent', figure1);
% 
% % Calculate limit - use max of multiple channels if plotting on same graph
% % Estimate x limit to try and avoid using too much CPU when drawing
% xlim(axes1, [0 (actualSampleInterval * maxSamples)]);
% 
% yRange = channelARangeMV + 0.5;
% ylim(axes1,[(-1 * yRange) yRange]);
% 
% hold(axes1,'on');
% grid(axes1, 'on');
% 
% title('PicoScope 4000 Series - Streaming Data Capture');
% xLabelStr = strcat('Time (', sampleIntervalTimeUnitsStr, ')');
% xlabel(xLabelStr);
% ylabel('Voltage (mV)');

% Get data values as long as power status has not changed (check for STOP button push inside loop)
i = 1;
while(hasAutoStopped == PicoConstants.FALSE && status.getStreamingLatestValues == PicoStatus.PICO_OK)
    
    ready = PicoConstants.FALSE;
   
    while(ready == PicoConstants.FALSE)

       status.getStreamingLatestValues = invoke(streamingGroupObj, 'getStreamingLatestValues');
       
       ready = invoke(streamingGroupObj, 'isReady');

       % Give option to abort from here
       flag = getappdata(gcf, 'run');
       drawnow;

       if(flag == 0)

            disp('STOP button clicked - aborting data collection.')
            break;

       end

       drawnow;

    end
    
    % Check for data
    [newSamples, startIndex] = invoke(streamingGroupObj, 'availableData');

    if (newSamples > 0)
        
        % Check if the scope has triggered
        [triggered, triggeredAt] = invoke(streamingGroupObj, 'isTriggerReady');

        if (triggered == PicoConstants.TRUE)

            % Adjust trigger position as MATLAB does not use zero-based
            % indexing
            
            bufferTriggerPosition = triggeredAt + 1;
            
            fprintf('Triggered - index in buffer: %d\n', bufferTriggerPosition);

            hasTriggered = triggered;

            % Adjust by 1 due to driver using zero indexing
            triggeredAtIndex = totalSamples + bufferTriggerPosition;

        end

        previousTotal = totalSamples;
        totalSamples = totalSamples + newSamples;

        % Printing to console can slow down acquisition - use for debug
%         fprintf('Collected %d samples, startIndex: %d total: %d.\n', newSamples, startIndex, totalSamples);
        
        % Position indices of data in buffer
        firstValuePosn = startIndex + 1;
        lastValuePosn = startIndex + newSamples;
        
        % Convert data values to milliVolts from the application buffers
%         bufferChAmV = adc2mv(pAppBufferChA.Value(firstValuePosn:lastValuePosn), channelARangeMV, maxADCCount);
%         bufferChBmV = adc2mv(pAppBufferChB.Value(firstValuePosn:lastValuePosn), channelBRangeMV, maxADCCount);
%         bufferChA = pAppBufferChA.Value(firstValuePosn:lastValuePosn);
%         bufferChB = pAppBufferChB.Value(firstValuePosn:lastValuePosn);
%         l = length(bufferChA);
%         bufferChAmV{i,1} = zeros(l, 1, 'int16');
%         bufferChBmV{i,1} = zeros(l, 1, 'int16');
            % adc2mv converts into double precision (64 bit)
            % data are in 16 bit...
%         bufferChAmV{i,1} = adc2mv(pAppBufferChA.Value(firstValuePosn:lastValuePosn), channelARangeMV, maxADCCount);
%         bufferChBmV{i,1} = adc2mv(pAppBufferChB.Value(firstValuePosn:lastValuePosn), channelBRangeMV, maxADCCount);
        bufferChAmV{i,1} = pAppBufferChA.Value(firstValuePosn:lastValuePosn);
        bufferChBmV{i,1} = pAppBufferChB.Value(firstValuePosn:lastValuePosn);
        i = i+1;
        % Process collected data further if required - this example plots
        % the data.
        
        % Copy data into final buffers
%         pBufferChAFinal.Value(previousTotal + 1:totalSamples) = bufferChAmV;
%         pBufferChBFinal.Value(previousTotal + 1:totalSamples) = bufferChBmV;
        
%         % Alberto
%         A(previousTotal + 1:totalSamples) = bufferChAmV;
%         B(previousTotal + 1:totalSamples) = bufferChBmV;
        
        % Time axis
        
        % Multiply by ratio mode as samples get reduced
%         time = (double(actualSampleInterval) * double(downSampleRatio)) * [previousTotal:(totalSamples - 1)];
        
        % Plot channel A only
%         plot(time, bufferChAmV);

        % Clear variables for use again
        clear bufferChAMV;
        clear firstValuePosn;
        clear lastValuePosn;
        clear startIndex;
        clear triggered;
        clear triggerAt;
   
    end
   
    % Check if auto stop has occurred
    hasAutoStopped = invoke(streamingGroupObj, 'autoStopped');

    if(hasAutoStopped == PicoConstants.TRUE)

       disp('AutoStop: TRUE - exiting loop.');
       break;

    end
   
    % Check if 'STOP' button pressed

    flag = getappdata(gcf, 'run');
%     drawnow;

    if(flag == 0)

        disp('STOP button clicked - aborting data collection.')
        break;
        
    end
 
end

% Close the STOP button window
if(exist('stopFig', 'var'))
    
    close('Stop Button');
    clear stopFig;
        
end

drawnow;

if(hasTriggered == PicoConstants.TRUE)
   
    fprintf('Triggered at overall index: %d\n', triggeredAtIndex);
    
end

% Take hold off the current figure
hold off;

fprintf('\n');

%% STOP THE DEVICE
% This function should be called regardless of whether auto stop is enabled
% or not.

status.stop = invoke(ps4000DeviceObj, 'ps4000Stop');

%% FIND THE NUMBER OF SAMPLES
% This is the number of samples held in the driver itself.
[status.noOfStreamingValues, numStreamingValues] = invoke(streamingGroupObj, 'ps4000NoOfStreamingValues');

fprintf('Number of samples available from the driver: %u.\n\n', numStreamingValues);

%% PROCESS DATA
% Process all data if required

% % Reduce size of arrays if required
% if(totalSamples < maxSamples)
%     
%     pBufferChAFinal.Value(totalSamples + 1:end) = [];
%     pBufferChBFinal.Value(totalSamples + 1:end) = [];
%  
% end

% channelAFinal = pBufferChAFinal.Value();
% channelBFinal = pBufferChBFinal.Value();

% % Plot total data on another figure
% finalFigure = figure;
% axes2 = axes('Parent', finalFigure);
% hold on;
% 
% title('Streaming Data Capture');
% xLabelStr = strcat('Time (', sampleIntervalTimeUnitsStr, ')');
% xlabel(xLabelStr);
% ylabel('Voltage (mV)');
% 
% % Find the maximum voltage range
% maxYRange = max(channelARangeMV, channelBRangeMV);
% ylim(axes2,[(-1 * maxYRange) maxYRange]);

% time = (double(actualSampleInterval) * double(downSampleRatio)) * [0:length(channelAFinal) - 1];
% plot(time, channelAFinal, time, channelBFinal);

% grid on;
% legend('Channel A', 'Channel B');


% hold off;
%% DISCONNECT DEVICE

disconnect(ps4000DeviceObj);



%% save
% Alberto

defAns = {datestr(now,'yyyymmdd')};
basename = inputdlg('save','Give basename for saving', 1, defAns);
filename = [basename{1,1}];

% check if file exists
cd(folder2save)
temp = dir('*.mat');
for i = 1:length(temp)
    if strcmp(filename, temp(i,1).name)
        warning('file already exits')
        warning('file not saved')
        return
    end
end

% save
save(filename, 'bufferChAmV', 'bufferChBmV', 'Tinterval', 'actualSampleInterval', 'channelSettings')
fprintf('\n\nSaved as %s\n', filename)

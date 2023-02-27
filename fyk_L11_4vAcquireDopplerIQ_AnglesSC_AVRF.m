
% Detect motion in ROI and aquire RF data
% Plane wave B mode for monitoring and detecting
% Moving Roi, DPIRoiOffset can be changed in GUI
% Two TPCs are used in this script
% Wait for triggering in before Doppler detecting
% Acquiring/processing/Saving automatically
% Changing the "EXPtimes" to set running turns


% Last update: Hu FOX, 22/06/2019,shenzhen for vantage 4.0



%% System parameters
clear;

filename = ('fyk_L11_4vAcquireDopplerIQ_AnglesSC');
pn = uigetdir;
Dnum = 1;
SavingNum = 10;

EXPtimes = 200;
ScaleMax = 500;  % scaling of the display function
ScaleMin = 0;  % scaling of the display function
ScaleRange = 0;  % scaling of the display function
ScaleOffset = 0;  % scaling of the display function
ne = 50;         % Set ne = number of detect acquisitions.
% DPIFrames   = 1;
BmodeFrames = 10;
DopplerFrames = 100;

maxVoltage = 61;
Fc = 6.25e6;
na  = 3;%number of views (angles) for compound imaging, can be even or odd
PRTus = 200; % in us
PRF = 1/(na*PRTus*1e-6);
dtheta = (12*pi/180)/(na);        %increment of angle:-2,0,2
startAngle = -(na-1)/2*dtheta;

WBFactor = 1.8;  %  control the element size for wide beam detect

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.simulateMode = 0;


% Specify Trans structure array.
Trans.name = 'L11-4v';
% Trans.frequency = 6.4286; % V1 nominal frequency in megahertz
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = maxVoltage;  % set maximum high voltage limit for pulser supply.

% wvl = Resource.Parameters.speedOfSound/(Trans.frequency*1e6);

%% Imaging parameters

% Specify Format structure array.
Format.transducer = 'L11-4v';   % 128 element linear array
Format.scanFormat = 'RLIN';     % rectangular linear array scan
Format.startDepth = 100;   % Acquisition depth in wavelengths150
Format.endDepth = 256;   % This should preferrably be a multiple of 128 samples.300

DPIFocusX = 0;
DPIFocusZ = Format.endDepth-(Format.endDepth-Format.startDepth)/2;%65;
DPIROI = [120 (Format.endDepth-Format.startDepth)/2];  % [width in element number (mush be even), depth in wavelength]

DPIRoiOffset.x = 0; % ROI offset in X direction
DPIRoiOffset.y = 0;

% Specify PData(1) structure array for Bmode Imaging
% PData(1).Format = 1;      % use first Format structure.
PData(1).PDelta = [Trans.spacing,0,0.5];
PData(1).Size(1) = ceil((Format.endDepth-Format.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((128*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*63.5,0,Format.startDepth]; % x,y,z of upper lft crnr.
PData(1).Region(1).Shape = struct(...
    'Name','Rectangle',...
    'Position',[0,0,Format.startDepth],...
    'width',ceil((128*Trans.spacing)/PData(1).PDelta(1)),...
    'height', Format.endDepth-Format.startDepth);
PData(1).Region(2).Shape = struct(...
    'Name','Rectangle',...
    'Position',[DPIFocusX,0,DPIFocusZ-DPIROI(2)/2],...
    'width', DPIROI(1),...
    'height', DPIROI(2));
PData(1).Region = computeRegions(PData(1));
% 
% % Specify PData(1) structure array for Shearwave visulization
% % PData(1).Format = 1;      % use first Format structure.
% PData(1).PDelta(1) = Trans.spacing;
% PData(1).PDelta(3) = 0.5;
% PData(1).Size(1) = ceil((Format.endDepth-Format.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
% PData(1).Size(2) = ceil((128*Trans.spacing)/PData(1).PDelta(1));
% PData(1).Size(3) = 1;      % single image page
% PData(1).Origin = [-Trans.spacing*63.5,0,Format.startDepth]; % x,y,z of upper lft crnr.
% PData(1).Region(1).Shape = struct(...
%     'Name','Rectangle',...
%     'Position',[DPIFocusX,0,DPIFocusZ-DPIROI(2)/2],...
%     'width', DPIROI(1),...
%     'height', DPIROI(2));
% PData(1).Region = computeRegions(PData(1));

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.function = 'movePoints';

maxAcqLength = sqrt(Format.endDepth^2 + (128*Trans.spacing)^2) - Format.startDepth;
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
% parameter of Rcvdata{2}
para.c = 1540;
para.fc = Fc;
para.fs = para.fc*4;
para.ts = 1/para.fs;
para.na = na;
para.ne = ne;
para.TXangle = startAngle+[0:na-1].*dtheta;
para.pitch = para.c/para.fc;
para.startDepth = Format.startDepth;
para.endDepth = Format.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128);
para.intervalSample = ceil(para.endDepth-para.startDepth)*4*2;
para.DetectTprf = PRTus*1e-6;
para.Tprf = PRTus*1e-6;   %?????????????xxxxxxxxx???????????????!!!!!!!!
para.prf = 1/para.Tprf;


% ROI Info
Depth = floor(PData(1).Region(2).Shape.height*2);
Width = floor(PData(1).Region(2).Shape.width);
Origin_x = floor(PData(1).Region(2).Shape.Position(1)-Width/2+65);
Origin_z = floor(((PData(1).Region(2).Shape.Position(3))-Format.startDepth)*2);
ROIinfo.Depth = Depth;
ROIinfo.Width = Width;
ROIinfo.Origin_x = Origin_x;
ROIinfo.Origin_z = Origin_z;

%% Specify Resources.
autofixed_sampleNum = 1024*ceil(para.intervalSample/1024);

% RcvBuffer for all raw data
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*autofixed_sampleNum;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;  % RcvBuffer is 64 cols using syn aper.
Resource.RcvBuffer(1).numFrames = BmodeFrames;     %

% InterBuffer for B Mode
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.InterBuffer(1).rowsPerFrame = PData(1).Size(1);
Resource.InterBuffer(1).colsPerFrame = PData(1).Size(2);

% ImageBuffer for reference Bmode image
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1); % this is for maximum depth
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);
Resource.ImageBuffer(1).numFrames = BmodeFrames;

% RcvBuffer for doppler wave raw data
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = ne*na*autofixed_sampleNum;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;  % RcvBuffer is 64 cols using syn aper.
Resource.RcvBuffer(2).numFrames = DopplerFrames;     %


% InterBuffer for DPI visualizaion
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;  % one intermediate buffer needed.
Resource.InterBuffer(2).rowsPerFrame = PData(1).Size(1);
Resource.InterBuffer(2).colsPerFrame = PData(1).Size(2);
Resource.InterBuffer(2).pagesPerFrame = ne;


Resource.DisplayWindow(1).Title = filename;
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');

DwWidth = ceil(PData(1).Region(1).Shape.width*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Region(1).Shape.height*2*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Colormap = gray(256);

%% Transmit parameters
% Specify Transmit waveform structure.
% - detect waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Fc/1e6,0.67,2,1];   % A, B, C, D

TW(2).type = 'parametric';
TW(2).Parameters = [Fc/1e6,0.67,8,1];   % A, B, C, D

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements),...
    'TXPD', [], ...
    'peakCutOff', 0.2, ...%for more signal change to 0.2 from 1.0
    'peakBLMax', 4.0), 1, na);



for ii = 1:na
    TX(ii).focus = 0.0;
    TX(ii).Steer = [(startAngle+(ii-1)*dtheta),0.0];%    TX(n).Steer = [(startAngle+(n-1)*dtheta),0.0];
    TX(ii).Apod = ones(1,Trans.numelements);
    TX(ii).Delay = computeTXDelays(TX(ii));
end


TPC(1).name = '2D';
TPC(1).maxHighVoltage = maxVoltage;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = 35;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structure,  be careful here
TGC.CntrlPts = [200,590,650,710,770,800,850,950];
TGC.rangeMax = Format.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);


%% Receive parameters
% Specify Receive structure arrays.
% - We need 2*na Receives for every frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = sqrt(Format.endDepth^2 + (128*Trans.spacing)^2) - Format.startDepth;
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
    'startDepth', Format.startDepth, ...
    'endDepth', Format.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'mode', 0, ...
    'callMediaFunc', 0), 1, na*BmodeFrames+ne*na*Resource.RcvBuffer(2).numFrames);

% - Set event specific Receive attributes for B-mode
for ii = 1:Resource.RcvBuffer(1).numFrames
    Receive(na*(ii-1)+1).callMediaFunc = 1;  % make media move per frame
    for j = 1:na  % numRays acquisitions per frame for BMode
        Receive(na*(ii-1)+j).Apod(1:128) = 1.0;
        Receive(na*(ii-1)+j).framenum = ii;
        Receive(na*(ii-1)+j).bufNum = 1;
        Receive(na*(ii-1)+j).acqNum = j;
    end

end

% - Set event specific Receive attributes for Doppler-mode
ind = na*BmodeFrames;
for N = 1:Resource.RcvBuffer(2).numFrames
    for j = 1:ne
        for k = 1:na
            Receive(ind+(j-1)*na+k).callMediaFunc = 1;
            Receive(ind+(j-1)*na+k).Apod(1:128) = 1.0;
            Receive(ind+(j-1)*na+k).framenum = N;
            Receive(ind+(j-1)*na+k).bufnum = 2;
            Receive(ind+(j-1)*na+k).acqNum = (j-1)*na+k;
        end
    end
    ind = ind + na*ne;
end



%% Reconstruction parameters
% Specify Recon structure arrays.

Recon = repmat(struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', zeros(1,1)), 1, 2);

Recon(1) = struct(...
    'senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'rcvBufFrame',-1, ...
    'IntBufDest', [1,1], ...  %inter buffer is only used for shearwave here
    'ImgBufDest', [1,-1], ...
    'RINums',(1:na)');


Recon(2) = struct(...
    'senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'rcvBufFrame',-1, ...
    'IntBufDest', [2,1], ...
    'ImgBufDest', [0,0], ...
    'RINums',(na+1:(ne*na+na))');


% Define ReconInfo structures.
% - ReconInfo for 2D frame.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'pagenum',1, ...
    'regionnum', 1), 1, na+ne*na);

ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
for j = 1:na
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
end
ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect

%  - ReconInfos for Doppler ensemble.
for j = 1:ne
    for ii = 1:na
        if ii == 1
            ReconInfo(na+(j-1)*na+ii).mode = 'replaceIQ';
        end
        ReconInfo(na+(j-1)*na+ii).txnum = ii;
        ReconInfo(na+(j-1)*na+ii).rcvnum = na*BmodeFrames+(j-1)*na+ii;  % ???????!!!!!!
        ReconInfo(na+(j-1)*na+ii).pagenum = j;
        ReconInfo(na+(j-1)*na+ii).regionnum = 2;
    end
end


%% Process parameters
% Specify Process structure array. (1) is used for B-mode imaging
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1, ...   % number of buffer to process.
    'framenum',-1, ...   % frame number in src buffer (-1 => lastFrame)
    'pdatanum',1, ...    % number of PData(1) structure (defines output figure).
    'norm',1, ...        % normalization method(1 means fixed)
    'pgain',3.0, ...            % pgain is image processing gain
    'persistMethod','simple', ...
    'persistLevel',20, ...
    'interp',1, ...      % method of interpolation (1=4pt interp)
    'compression',0.5, ...      % X^0.5 normalized to output word size
    'reject',2,...
    'mappingMode','full', ...
    'display',1, ...     % display image after processing
    'displayWindow',1};


Process(2).classname = 'External';
Process(2).method = 'SavingIQData';
Process(2).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',2,...
    'srcframenum',1,...
    'dstbuffer','none'};
    
%% SeqControl and Events
% - Noop to allow time for charging external cap.
SeqControl(1).command = 'noop';
SeqControl(1).argument = 500000; % wait 100 msec.

% - time between detect acquisitions
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = PRTus;               % 100 usec (10KHz),500 usec (2 KHz)
TTNAQ = 2;

% - Return to Matlab
SeqControl(3).command = 'returnToMatlab';
RTML = 3;

SeqControl(4).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions     
SeqControl(4).argument = Format.endDepth*2/(para.c)+900;
TTNLine = 4;

SeqControl(5).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(5).argument = 30000;  % 30000 usec = 30msec time between frames
TTNF = 5;

% -- Change to Profile 2 (Doppler)
SeqControl(6).command = 'setTPCProfile';
SeqControl(6).condition = 'next';
SeqControl(6).argument = 2;
TCPDP = 6;

% -- Change to Profile 1 (2D)
SeqControl(7).command = 'setTPCProfile';
SeqControl(7).condition = 'next';
SeqControl(7).argument = 1;
TCPBB = 7;

SeqControl(8).command = 'triggerOut';
SeqControl(8).condition = 'syncADC_CLK'; % input BNC #1 falling edge
TGOUT = 8;

nsc = 9;

% Specify Event structure arrays.
n = 1;

nStartBMode = n;

for ii = 1:Resource.RcvBuffer(1).numFrames   %Resource.RcvBuffer(1).numFrames = BmodeFrames=10;
    
    for j = 1: na
        Event(n).info = 'B Mode';
        Event(n).tx = j;  %PWNum*PWFocusXNum+j;         % use 1st TX structure.PWNum = 1
        Event(n).rcv = na*(ii-1)+j;      % use 2nd Rcv structure.numRays = 48;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = [TTNLine];% time between frames, SeqControl struct defined below.
        n = n+1;
        
    end
    Event(n-1).seqControl = TTNF;
 
    Event(n).info = 'transferToHost';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon =0;      % reconstruction
    Event(n).process = 0;    % processing
    Event(n).seqControl = nsc; 
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n + 1;
    
    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % processing
    Event(n).seqControl = 0; % return to Matlab
    n = n+1;
    
    
    Event(n).info = 'ext func to draw ROI';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 0;    % process
    Event(n).seqControl = RTML;
    n = n+1;
    
end


Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = nsc; % jump command
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStartBMode;
nsc = nsc + 1;
n = n+1;

nStartDoppler = n;  
%

Event(n).info = 'TPC Doppler';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = TCPDP;
n = n+1;

ind = na*BmodeFrames;
for N = 1:Resource.RcvBuffer(2).numFrames
    for j = 1:ne
        for ii=1:na
            Event(n).info = 'Acquire data';
            Event(n).tx = ii; %  PWNum*PWFocusXNum+numRays+1;
            Event(n).rcv = ind+na*(j-1)+ii; %numRays*BmodeFrames+j;%%%
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = [TTNAQ]; %100us per detect
            n = n+1;
        end
    end
    Event(n-1).seqControl = [TTNF,nsc,RTML];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;


    Event(n).info = 'reconstruct to IQ';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 2;      % reconstruction
    Event(n).process = 0;    % no process
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Save one frame data';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no recon
    Event(n).process = 4;    % save
    Event(n).seqControl = 0;
    n = n+1;

    ind = ind + na*ne;
end


Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = nsc; % jump command
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStartDoppler+1;
nsc = nsc + 1;
n = n+1;



%% User specified UI Control Elements

% Define UIPos, which contains the default GUI positions - three columns of 10 controls. The x,y
%    locations increment up columns, with each column being a separate page. The origin
%    specified by UIPos is the lower left corner of a virtual box that encloses the control.
UIPos = zeros(10,2,3);
UIPos(:,1,1) = 0.0625;
UIPos(:,1,2) = 0.375;
UIPos(:,1,3) = 0.6875;
UIPos(:,2,1) = 0.0:0.1:0.9;
UIPos(:,2,2) = 0.0:0.1:0.9;
UIPos(:,2,3) = 0.0:0.1:0.9;

% Define slider group offsets and sizes. All units are normalized.
SG = struct('TO',[0.0,0.0975],...   % title offset
    'TS',[0.25,0.025],...   % title size
    'TF',0.8,...            % title font size
    'SO',[0.0,0.06],...     % slider offset
    'SS',[0.25,0.031],...   % slider size
    'EO',[0.075,0.031],...   % edit box offset
    'ES',[0.11,0.031]);     % edit box size




% The follow UIs are not using User##
Pos = UIPos(2,:,3);



% Measure
UI(14).Control = {'Style','pushbutton',...
    'String','Measure',...
    'Units','normalized',...
    'Position',[Pos+[-0.13 0.05],0.12,0.04],...
    'FontUnits','normalized',...
    'FontSize',0.5,...
    'Callback',{@measure}};
UI(14).Callback = text2cell('%-UI#14Callback');

UI(15).Control = {'UserC3','Style','VsSlider','Label','ScaleRange',...
    'SliderMinMaxVal',[-50,50,ScaleRange],...
    'SliderStep',[0.01,0.01],'ValueFormat','%3.0f'};
UI(15).Callback = text2cell('%-UI#15Callback');

UI(16).Control = {'UserC2','Style','VsSlider','Label','ScaleOffset',...
    'SliderMinMaxVal',[-50,50,ScaleOffset],...
    'SliderStep',[0.01,0.01],'ValueFormat','%3.0f'};
UI(16).Callback = text2cell('%-UI#16Callback');
% 
VPAdj = 1;
UI(17).Control = {'UserC4','Style','VsButtonGroup','Title','Velocity/Power Adjustment',...
    'NumButtons',2,'Labels',{'Velocity','Power'}};
UI(17).Callback = {'assignin(''base'',''VPAdj'',UIState)'};


% External function definitions

% Create the External Processing function, using the new function handle
% approach. Note that in this example no import was used. The user can
% directly call the new approach functions using the full package path
% You can use the TAB-key after you place a '.' to get suggestions for
% available functions. 
%
% Note the callback function is a real function handle as indicated by the
% @ sign. You can now mark this function handle, and right click and then
% select 'Open runTimeMon' and the editor will jump to the function
% definition in your script
EF(1).Function = vsv.seq.function.ExFunctionDef('savingIQData', @()savingIQData(obj1,obj2));



% Handle for Shearwave figure
DPIHandle = figure('Name','ShearWaveVisulization',...
    'NumberTitle','off','Visible','off',...
    'Position',[Resource.DisplayWindow(1).Position(1)+250, ... % left edge
    Resource.DisplayWindow(1).Position(2), ... % bottom
    [DPIROI(1)+3,DPIROI(2)]*7], ...            % width, height
    'CloseRequestFcn',{@closeIQfig});


% figClose is used for figclose function in shearwave visualization
figClose = 0;

% adjStatus is used for checking ROI or focus adjustment
adjStatus = 1;

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = ne;


% Save all the structures to a .mat file.
save(['MatFiles\',filename]);

return

%% **** Callback routines to be converted by text2cell function. ****

%-UI#11Callback - ROI Width change
TX = evalin('base','TX');
UI = evalin('base','UI');
PData = evalin('base','PData');
Trans = evalin('base','Trans');
DPIROI = evalin('base','DPIROI');
Format = evalin('base','Format');
DPIHandle = evalin('base','DPIHandle');
na = evalin('base','numRays');
WBFactor = evalin('base','WBFactor');

DPIFocusX = evalin('base','DPIFocusX');
DPIFocusZ = evalin('base','DPIFocusZ');


if mod(UIValue,2)==1
    DPIROI(1) = UIValue+1;
    set(UI(11).handle(3),'String',UIValue+1);
else
    DPIROI(1) = UIValue;
end
assignin('base','DPIROI',DPIROI);

x = DPIFocusX;
z = DPIFocusZ;

% change PData(1) size
PData(1).Region(2).Shape = struct(...
    'Name','Rectangle',...
    'Position',[DPIFocusX,0,DPIFocusZ-DPIROI(2)/2],...
    'width', DPIROI(1),...
    'height', DPIROI(2));
PData(1).Region = computeRegions(PData(1));
assignin('base','TX',TX);
assignin('base','PData',PData);


Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX','PData'};
assignin('base','Control', Control);
assignin('base','adjStatus',1);

% pos = get(DPIHandle,'Position')
% set(DPIHandle,'Position',[pos(1:2),[DPIROI(1)+3,DPIROI(2)]*7]);


return
%-UI#11Callback


%-UI#12Callback - ROI Hight change
TX = evalin('base','TX');
PData = evalin('base','PData');
Trans = evalin('base','Trans');
DPIROI = evalin('base','DPIROI');
Format = evalin('base','Format');
DPIHandle = evalin('base','DPIHandle');
na = evalin('base','numRays');
WBFactor = evalin('base','WBFactor');

DPIFocusX = evalin('base','DPIFocusX');
DPIFocusZ = evalin('base','DPIFocusZ');


DPIROI(2) = UIValue;
assignin('base','DPIROI',DPIROI);

x = DPIFocusX;
z = DPIFocusZ;

% change PData(1) size
PData(1).Region(2).Shape = struct(...
    'Name','Rectangle',...
    'Position',[DPIFocusX,0,DPIFocusZ-DPIROI(2)/2],...
    'width', DPIROI(1),...
    'height', DPIROI(2));
PData(1).Region = computeRegions(PData(1));
assignin('base','TX',TX);
assignin('base','PData',PData);


Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Format','TX','PData','Receive','Recon','DisplayWindow'};
assignin('base','Control', Control);
assignin('base','adjStatus',1);

pos = get(DPIHandle,'Position')
set(DPIHandle,'Position',[pos(1:2),[DPIROI(1)+3,DPIROI(2)]*7]);


return
%-UI#12Callback


%-UI#14Callback
measure(varargin)

UI = evalin('base','UI');
freeze = evalin('base','freeze');
figClose = evalin('base','figClose');

% if freeze == 0 && figClose == 0
if figClose == 0

    % change the start event
%     CurrentState = get(UI(14).handle,'String')
% 
%     if strcmp(CurrentState,'Measure')
        
        nStart = evalin('base','nStartDoppler');
        Control = evalin('base','Control');
        Control.Command = 'set&Run';
        Control.Parameters = {'Parameters',1,'startEvent',nStart};
        evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
        assignin('base','Control',Control);
        set(UI(14).handle,'String','Monitor')
        
%     else
% 
%         nStart = 1;
%         Control = evalin('base','Control');
%         Control.Command = 'set&Run';
%         Control.Parameters = {'Parameters',1,'startEvent',nStart};
%         evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
%         assignin('base','Control',Control);
%         set(UI(14).handle,'String','Measure')
%         
%                 disp('yes')
%     end
    
end
return
%-UI#14Callback

%-UI#15Callback - ScaleRange change
ScaleMax = evalin('base','ScaleMax');
ScaleMin = evalin('base','ScaleMin');
ScaleRange = evalin('base','ScaleRange');
ScaleOffset = evalin('base','ScaleOffset');

ScaleRange = UIValue;
ScaleMax = sign(ScaleOffset)*2^abs(ScaleOffset) + 2^ScaleRange;
ScaleMin = sign(ScaleOffset)*2^abs(ScaleOffset) - 2^ScaleRange;
disp([' ScaleMin = ',num2str(ScaleMin),'ScaleMax = ',num2str(ScaleMax)]);

assignin('base','ScaleMax',ScaleMax);
assignin('base','ScaleMin',ScaleMin);
assignin('base','ScaleRange',ScaleRange);
assignin('base','ScaleOffset',ScaleOffset);

return

%-UI#16Callback - ScaleOffset change
ScaleMax = evalin('base','ScaleMax');
ScaleRange = evalin('base','ScaleRange');
ScaleOffset = evalin('base','ScaleOffset');

ScaleOffset = UIValue;
ScaleMax = sign(ScaleOffset)*2^abs(ScaleOffset) + 2^ScaleRange;
ScaleMin = sign(ScaleOffset)*2^abs(ScaleOffset) - 2^ScaleRange;
disp([' ScaleMin = ',num2str(ScaleMin),'ScaleMax = ',num2str(ScaleMax)]);

assignin('base','ScaleMax',ScaleMax);
assignin('base','ScaleMin',ScaleMin);
assignin('base','ScaleRange',ScaleRange);
assignin('base','ScaleOffset',ScaleOffset);

return

%% External functions


function SavingIQData(IBuffer,QBuffer)
%
    Dnum = evalin('base','Dnum');
    para = evalin('base','para');
    ROIinfo = evalin('base','ROIinfo');
    pn = evalin('base','pn');
    SavingNum = evalin('base','SavingNum');
    
    if Dnum <= SavingNum
        if ~isempty(pn) % fn will be zero if user hits cancel
            savefast([pn,'\Data',int2str(Dnum),'.mat'],'IBuffer','QBuffer','para','ROIinfo');
            fprintf('The %dth data has been saved!\n ',Dnum);
            Dnum = Dnum+1;
            assignin('base','Dnum',Dnum);
        else
            disp('The data is not saved.');
        end
    end

end


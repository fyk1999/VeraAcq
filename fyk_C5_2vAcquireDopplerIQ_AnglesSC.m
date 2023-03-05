

%% System parameters
clear;

filename = ('fyk_C5_2vAcquireDopplerIQ_AnglesSC');
pn = 'E:\FYK\Delete';


ScaleMax = 500;  % scaling of the display function
ScaleMin = 0;  % scaling of the display function
ScaleRange = 0;  % scaling of the display function
ScaleOffset = 0;  % scaling of the display function
ne = 50;         % Set ne = number of detect acquisitions.
% DPIFrames   = 1;
BmodeFrames = 10;
DopplerFrames = 100;

maxVoltage = 50;
Fc = 6.25e6;
na  = 3;%number of views (angles) for compound imaging, can be even or odd
numAccum = 2;
PRFdoppler = 1e3;
PRF = PRFdoppler*na*numAccum;
PRT = 1/PRF;
dtheta = (12*pi/180)/(na);        %increment of angle:-2,0,2
startAngle = -(na-1)/2*dtheta;


% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.simulateMode = 1;


% Specify Trans structure array.
Trans.name = 'C5-2v';
% Trans.frequency = 6.4286; % V1 nominal frequency in megahertz
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = maxVoltage;  % set maximum high voltage limit for pulser supply.
radius = Trans.radius;
scanangle = Trans.numelements*Trans.spacing/radius;
theta = -(scanangle/2); % angle to left edge from centerline


%% Imaging parameters

% Specify Format structure array.
Format.transducer = 'L11-4v';   % 128 element linear array
Format.scanFormat = 'RLIN';     % rectangular linear array scan
Format.startDepth = 20;   % Acquisition depth in wavelengths150
Format.endDepth = 256;   % This should preferrably be a multiple of 128 samples.300
P.startDepth = Format.startDepth;
P.endDepth = Format.endDepth;



DPIROI.focus = Format.startDepth + (Format.endDepth-Format.startDepth)/2;  % center point at default;
DPIROI.angles = scanangle/2;
DPIROI.height  = (Format.endDepth-Format.startDepth)/2;   


% Specify PData structure array.
PData(1).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData(1).Origin(1,2) = 0;
PData(1).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;
PData(1).Region.Shape = struct(...
                   'Name','Sector',...
                   'Position',[0,0,-radius],...
                   'r1',radius+P.startDepth,...
                   'r2',radius+P.endDepth,...
                   'angle',scanangle);
PData(1).Region(2).Shape = struct(...
                   'Name','Sector',...
                   'Position',[0,0,-radius],...
                   'r1',radius+DPIROI.focus - DPIROI.height/2,...
                   'r2',radius+DPIROI.focus + DPIROI.height/2,...
                   'angle',DPIROI.angles);
PData(1).Region = computeRegions(PData(1));




% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.MP(1,:) = [0,0,10,1.0];
Media.MP(2,:) = [(radius+10)*sin(-0.2608),0,(radius+10)*cos(-0.2608)-radius,1.0];
Media.MP(3,:) = [(radius+10)*sin(0.2608),0,(radius+10)*cos(0.2608)-radius,1.0];
Media.MP(4,:) = [(radius+10)*sin(-0.5267),0,(radius+10)*cos(-0.5267)-radius,1.0];
Media.MP(5,:) = [(radius+10)*sin(0.5267),0,(radius+10)*cos(0.5267)-radius,1.0];
Media.MP(6,:) = [0,0,40,1.0];
Media.MP(7,:) = [0,0,70,1.0];
Media.MP(8,:) = [(radius+70)*sin(-0.2608),0,(radius+70)*cos(-0.2608)-radius,1.0];
Media.MP(9,:) = [(radius+70)*sin(0.2608),0,(radius+70)*cos(0.2608)-radius,1.0];
Media.MP(10,:) = [(radius+70)*sin(-0.5267),0,(radius+70)*cos(-0.5267)-radius,1.0];
Media.MP(11,:) = [(radius+70)*sin(0.5267),0,(radius+70)*cos(0.5267)-radius,1.0];
Media.MP(12,:) = [0,0,100,1.0];
Media.MP(13,:) = [0,0,130,1.0];
Media.MP(14,:) = [(radius+130)*sin(-0.2608),0,(radius+130)*cos(-0.2608)-radius,1.0];
Media.MP(15,:) = [(radius+130)*sin(0.2608),0,(radius+130)*cos(0.2608)-radius,1.0];
Media.MP(16,:) = [(radius+130)*sin(-0.5267),0,(radius+130)*cos(-0.5267)-radius,1.0];
Media.MP(17,:) = [(radius+130)*sin(0.5267),0,(radius+130)*cos(0.5267)-radius,1.0];
Media.MP(18,:) = [0,0,160,1.0];
Media.MP(19,:) = [0,0,190,1.0];
Media.MP(20,:) = [(radius+190)*sin(-0.2608),0,(radius+190)*cos(-0.2608)-radius,1.0];
Media.MP(21,:) = [(radius+190)*sin(0.2608),0,(radius+190)*cos(0.2608)-radius,1.0];
Media.MP(22,:) = [(radius+190)*sin(-0.5267),0,(radius+190)*cos(-0.5267)-radius,1.0];
Media.MP(23,:) = [(radius+190)*sin(0.5267),0,(radius+190)*cos(0.5267)-radius,1.0];
Media.function = 'movePoints';



%% The Parameters that need to be saved
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
para.DetectTprf = 0;
para.Tprf = 0;   %?????????????xxxxxxxxx???????????????!!!!!!!!
para.prf = PRF;


% ROI Info
ROIinfo = ["Unable to calculate for now!"];



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


% DisplayWindow for B-mode image show
Resource.DisplayWindow(1).Title = filename;
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Matlab';
%Resource.DisplayWindow(1).numFrames = 20;
%Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);
Resource.VDAS.dmaTimeout = 2000;

% Resource.DisplayWindow(1).Title = filename;
% Resource.DisplayWindow(1).pdelta = 0.35;
% ScrnSize = get(0,'ScreenSize');
% DwWidth = ceil(PData(1).Region(1).Shape.width*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
% DwHeight = ceil(PData(1).Region(1).Shape.height*2*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
% Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
%     DwWidth, DwHeight];
% Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
% Resource.DisplayWindow(1).Colormap = gray(256);
% Resource.DisplayWindow(1).Type = 'Matlab';
% Resource.VDAS.dmaTimeout = 2000;





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
    'callMediaFunc', 0), 1, na*BmodeFrames+ne*na*Resource.RcvBuffer(2).numFrames*2);

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
            % first acquistion
            Receive(ind+(j-1)*na*2+(k-1)*2+1).callMediaFunc = 1;
            Receive(ind+(j-1)*na*2+(k-1)*2+1).Apod(1:128) = 1.0;
            Receive(ind+(j-1)*na*2+(k-1)*2+1).framenum = N;
            Receive(ind+(j-1)*na*2+(k-1)*2+1).bufnum = 2;
            Receive(ind+(j-1)*na*2+(k-1)*2+1).acqNum = (j-1)*na+k;
            Receive(ind+(j-1)*na*2+(k-1)*2+1).mode = 0;
            % accumulate acquistion
            Receive(ind+(j-1)*na*2+(k-1)*2+2).callMediaFunc = 0;
            Receive(ind+(j-1)*na*2+(k-1)*2+2).Apod(1:128) = 1.0;
            Receive(ind+(j-1)*na*2+(k-1)*2+2).framenum = N;
            Receive(ind+(j-1)*na*2+(k-1)*2+2).bufnum = 2;
            Receive(ind+(j-1)*na*2+(k-1)*2+2).acqNum = (j-1)*na+k;
            Receive(ind+(j-1)*na*2+(k-1)*2+2).mode = 1;
        end
    end
    ind = ind + na*ne*2;
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
pers = 20;
cmpFactor = 40;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',1.0,...     % pgain is image processing gain
    'reject',2,...
    'grainRemoval','none',...
    'persistMethod','none',...
    'persistLevel',pers,...
    'interpMethod','4pt',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',cmpFactor,...
    'mappingMethod','full',...
    'display',1,...      % display image after processing
    'displayWindow',1};


Process(2).classname = 'External';
Process(2).method = 'saveIQData';
Process(2).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',2,...
    'srcframenum',1,...
    'dstbuffer','none'};


Process(3).classname = 'External';
Process(3).method = 'drawROI';
Process(3).Parameters = {'srcbuffer','none'};


Process(4).classname = 'External';
Process(4).method = 'debugFunc';
Process(4).Parameters = {'srcbuffer','none'};



%% SeqControl and Events
% - Noop to allow time for charging external cap.
SeqControl(1).command = 'noop';
SeqControl(1).argument = 500000; % wait 100 msec.

% - time between detect acquisitions
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = round(PRT*1e6);               % 100 usec (10KHz),500 usec (2 KHz)
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
TCPD = 6;

% -- Change to Profile 1 (2D)
SeqControl(7).command = 'setTPCProfile';
SeqControl(7).condition = 'next';
SeqControl(7).argument = 1;
TCPBB = 7;


SeqControl(8).command = 'loopCnt';
SeqControl(8).argument = numAccum;
SLCN = 8;


SeqControl(9).command = 'sync';
SYNC = 9;


nsc = 10;

% Specify Event structure arrays.
n = 1;

nStartMonitor = n;

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
    Event(n).process = 3;    % process
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
SeqControl(nsc).argument = nStartMonitor;
nsc = nsc + 1;
n = n+1;


%%%%%% Another Loop %%%%%%

nStartMeasure = n; 
%
Event(n).info = 'TPC Doppler';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = TCPD;
n = n+1;


Event(n).info = 'sync after branch';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = SYNC;
n = n+1;


ind = na*BmodeFrames;
for N = 1:Resource.RcvBuffer(2).numFrames
    for j = 1:ne
        for ii=1:na
            for k = 1:numAccum
                if(k==1)
                    Event(n).info = '1st acquisition';
                    Event(n).tx = ii; %  PWNum*PWFocusXNum+numRays+1;
                    Event(n).rcv = ind+na*(j-1)*2+(ii-1)*2+1; %numRays*BmodeFrames+j;%%%
                    Event(n).recon = 0;      % no reconstruction.
                    Event(n).process = 0;    % no processing
                    Event(n).seqControl = [TTNAQ]; %100us per detect
                    n = n+1;
                else
                    Event(n).info = 'Accumulate acquisition';
                    Event(n).tx = ii; %  PWNum*PWFocusXNum+numRays+1;
                    Event(n).rcv = ind+na*(j-1)*2+(ii-1)*2+2; %numRays*BmodeFrames+j;%%%
                    Event(n).recon = 0;      % no reconstruction.
                    Event(n).process = 0;    % no processing
                    Event(n).seqControl = [TTNAQ]; %100us per detect
                    n = n+1;
                end
            end
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
    Event(n).process = 2;    % save
    Event(n).seqControl = 0;
    n = n+1;

    ind = ind + na*ne*2;
end


Event(n).info = 'Jump back to first event of Measure';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = nsc; % jump command
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStartMeasure+1;
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

import vsv.seq.uicontrol.*;

% - cBranch Gui Button
UI(1).Control = VsButtonControl('LocationCode','UserB2',...
    'Label','Measure','Callback',@measure);


% ROI adjustment
UI(2).Control = VsSliderControl('LocationCode','UserB5', ...
    'Label','ROI angle',...
    'SliderMinMaxVal',[0.4,0.8,DPIROI.angles],...
    'SliderStep',[2/80,10/80],'ValueFormat','%3.0f',...
    'Callback',@angleChange);


UI(3).Control = VsSliderControl('LocationCode','UserB4', ...
    'Style','VsSlider','Label','ROI Hight',...
    'SliderMinMaxVal',[20,(Format.endDepth-Format.startDepth),DPIROI.height],...
    'SliderStep',[2/100,10/100],'ValueFormat','%3.0f',...
    'Callback',@heightChange);


% Focus adjustment
checkFocusAdj = 1;
UI(4).Control = vsv.seq.uicontrol.VsButtonGroupControl('LocationCode','UserB3',...
    'Title','Adjust Focus',...
    'NumButtons',2,...
    'PossibleCases',   {'off','on'},...
    'Callback', @AdjustFocusCallback);



%% User specified External function

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
EF(1).Function = vsv.seq.function.ExFunctionDef('saveIQData', @saveIQData);
EF(2).Function = vsv.seq.function.ExFunctionDef('drawROI', @drawROI);
EF(3).Function = vsv.seq.function.ExFunctionDef('debugFunc', @debugFunc);

%% other thing

% Handle for Doppler figure
%dpiHandle = figure('Name','DopplerModeVisulization',...
%     'NumberTitle','off','Visible','off',...
%     'Position',[Resource.DisplayWindow(1).Position(1)+250, ... % left edge
%     Resource.DisplayWindow(1).Position(2), ... % bottom
%     [DPIROI(1)+3,DPIROI(2)]*7], ...            % width, height
%     'CloseRequestFcn',{@closeIQfig});

bmodeHandle = [];

% figClose is used for figclose function in shearwave visualization
figClose = 0;

% adjStatus is used for checking ROI or focus adjustment
adjStatus = 1;

%% Save .Mat file

% Save all the structures to a .mat file.
save(['MatFiles\',filename]);
VSX
return



%% External functions

% - save one frame of iq buffer data
function saveIQData(IData,QData)
%
    para = evalin('base','para');
    DPIROI = evalin('base','DPIROI');
    pn = evalin('base','pn');

    persistent SavingNum;
    persistent Dnum ;
    if isempty(Dnum)
        Dnum = 1;
    end
    if isempty(SavingNum)
        SavingNum = 5;           %the frame number that will be saved
    end

    if Dnum <= SavingNum
        if ~isempty(pn) % fn will be zero if user hits cancel
            savefast([pn,'\Data',int2str(Dnum),'.mat'],'IData','QData','para','DPIROI');
            fprintf('The %dth data has been saved!\n ',Dnum);
            Dnum = Dnum+1;
            assignin('base','Dnum',Dnum);
        else
            disp('The data is not saved.');
        end
    end

end



% Draw ROI of doppler in the image
function drawROI()

global recHandle1 markHandle

figClose  = evalin('base','figClose');
bmodeHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
DPIROI = evalin('base','DPIROI');
radius = evalin('base','radius');

switch figClose    
    case 1  % not freeze status, but fig has been closed
        assignin('base','figClose',2);
        delete(gcf)
        return
    case 2  % fig was closed at freeze status
        return
    case 0         
        adjStatus = evalin('base','adjStatus');            
        if adjStatus == 1        
            % Mark on bmode figure
            if ishandle(bmodeHandle)
                if  isempty(recHandle1) || ~ishandle(recHandle1)
                    figure(bmodeHandle), hold on,
                    %recHandle1 = rectangle('Position',[x-DPIROI(1)/2,z-DPIROI(2)/2,DPIROI(1),DPIROI(2)],'EdgeColor','r');
                    [newX , newY] = drawSector(DPIROI,radius);
                    recHandle1 = plot(newX,newY,'r-');  %plot is called for only once
                    markHandle = plot(0,DPIROI.focus,'xr','MarkerFaceColor','r','MarkerSize',8,'Linewidth',2);hold off;
                else
                    [newX , newY] = drawSector(DPIROI,radius);
                    set(recHandle1,'XData',newX,'YData',newY);
                    set(markHandle,'YData',DPIROI.focus);
                end
            end                       
        end            
end
assignin('base','BmodeHandle',bmodeHandle);

end



% - debug external function
function debugFunc()
    bmodeHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');    
    disp(bmodeHandle);
    pause
end


%% Callback functions

% - branch from one event loop to another loop
function measure(~, ~)
    nStart = evalin('base','nStartMeasure');
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Parameters',1,'startEvent',nStart};
    %evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% - Close handle
function closeIQfig(varargin)

exitValue = evalin('base','vsExit');
freeze = evalin('base','freeze');
% if GUI is not freeze, just assign value 1 for EF to close the fig
assignin('base','figClose',1);
if exitValue  % if GUI has been closed, delete fig
    delete(gcf)
end
if freeze  % if GUI is freeze, delete and assign 2 to figClose
    assignin('base','figClose',2);
    delete(gcf);
end

end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% - change ROI width
function angleChange(~,~,UIValue)
PData = evalin('base','PData');
DPIROI = evalin('base','DPIROI');
radius = evalin('base','radius');
if mod(UIValue,2)==1    % must be even
    DPIROI.angles = UIValue+1;
else
    DPIROI.angles = UIValue;
end
assignin('base','DPIROI',DPIROI);

% change PData(1).region(2).shape
PData(1).Region(2).Shape.angle = DPIROI.angles;
PData(1).Region = computeRegions(PData(1));
assignin('base','PData',PData);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData'};
assignin('base','Control', Control);
assignin('base','adjStatus',1);

end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% - ROI height change
function heightChange(~,~,UIValue)
PData = evalin('base','PData');
DPIROI = evalin('base','DPIROI');
radius = evalin('base','radius');

DPIROI.height = UIValue;
assignin('base','DPIROI',DPIROI);

% change PData(1) size
PData(1).Region(2).Shape.r1 = radius + DPIROI.focus - DPIROI.height/2;
PData(1).Region(2).Shape.r2 = radius + DPIROI.focus + DPIROI.height/2;
PData(1).Region = computeRegions(PData(1));
assignin('base','PData',PData);
%Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData'};
assignin('base','Control', Control);
assignin('base','adjStatus',1);
%pos = get(dpiHandle,'Position');
%set(dpiHandle,'Position',[pos(1:2),[DPIROI(1)+3,DPIROI(2)]*7]);

end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -  Unkown function
function focusAdj(varargin)

freeze = evalin('base','freeze');
checkFocusAdj = evalin('base','checkFocusAdj');
bmodeHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
bmodeAxes = get(bmodeHandle,'currentAxes');
DPIROI = evalin('base','DPIROI');
radius = evalin('base','radius');
Fc = evalin('base','Fc');

if checkFocusAdj == 2 && freeze == 0
    pos = get(bmodeAxes,'CurrentPoint');
    DPIROI.focus = round(pos(3));
    assignin('base','adjStatus',1);
    PData = evalin('base','PData');
    PData(1).Region(2).Shape.r1 = radius + DPIROI.focus - DPIROI.height/2;
    PData(1).Region(2).Shape.r2 = radius + DPIROI.focus + DPIROI.height/2;
    PData(1).Region = computeRegions(PData(1));
    assignin('base','PData',PData);
    assignin('base','DPIROI',DPIROI);
    %Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData'};
    assignin('base','Control', Control);
    fprintf("Focus Z = %0.2f(cm)\n",100*DPIROI.focus*1540/Fc);
end

end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% - 
function AdjustFocusCallback(~,~,UIValue)

Resource = evalin('base','Resource');
% checkFocusAdj =  evalin('base','checkFocusAdj');
checkFocusAdj = UIValue;
fh = Resource.DisplayWindow(1).figureHandle;
% 'set(Resource.DisplayWindow(1).figureHandle,''WindowButtonDownFcn'',{@focusAdj,''focusAdj''});';

if UIValue == 2  
    set(fh,'WindowButtonDownFcn',{@focusAdj,'focusAdj'});
else
    set(fh,'WindowButtonDownFcn',[]);
end

assignin('base','checkFocusAdj',checkFocusAdj);

end



%% Built-in function
function [x,y]  = drawSector(ROI,radius)  
x0 = 0;
y0 = -radius;
theta = ROI.angles;
a1 = 0.5*pi-0.5*theta;
a2 = a1+theta;
r1 = ROI.focus-ROI.height/2 - y0;
r2 = ROI.focus+ROI.height/2 - y0;
t = linspace(a1,a2);
x = x0 + r1*cos(t);
y = y0 + r1*sin(t);
t = flip(t);
x2 = x0 + r2*cos(t);
y2 = y0 + r2*sin(t);
x = [x x2 x(1)];
y = [y y2 y(1)];
end

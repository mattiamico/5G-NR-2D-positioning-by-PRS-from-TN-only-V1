close all ; clear all


formatOut = 'yyyy-mm-dd_HH.MM';
timestamp = datestr(now,formatOut);
logFileName = strcat('logFile_',timestamp,'.txt');
diary(logFileName) 

%% FIXME?: create output files folder with name: "<timestamp>-<rootName>" 
outputDirectory = strcat('./',timestamp);
mkdir(outputDirectory);   % create the directory

savePlots = true


%%%%%% Simulation Parameters %%%%%%

nFrames = 1; % Number of 10 ms 5G-NR frames ( 1 FRAME = 10 SUBFRAMES ) 
fc = 3e9; % Carrier frequency [Hz] ( ie 3 [GHz] )


validationMode = true ; % No path loss attenuations from propagation. No TOA derivation from correlation


%% Specify UE and gNB positions

% Configure UE position in xy-coordinate plane
UEPos = [500 -20]; % [m]
%UEPos = [17000 17000]; % [m]

% Configure number of gNBs and locate them at random positions in xy-coordinate plane
numgNBs = 4;
rng('default');  % Set RNG state for repeatability

%gNBPos = getgNBPositions(numgNBs); % [m] (ORIGINAL MATLAB TUTORIAL) 
r1 = 4*10^3; % [m]
r2 = 6*10^3; % [m] % == r1 for validation

angularlyEquispaced = false

%central_position = [0,0];
central_position = UEPos;

gNBPos = getgNB2DPositions(numgNBs, r1, r2, central_position, angularlyEquispaced);  % [m]


disp('gNBs positions [m]')
celldisp(gNBPos)


% Plot UE and gNB positions
%plotgNBAndUE_2DPositions(gNBPos,UEPos,1:numgNBs,r1,r2);
plotgNBAndUE_2DPositions(gNBPos,UEPos,1:numgNBs,r1,r2,central_position, 'TRUE Horizontal positions', savePlots,outputDirectory);


%%%%%% Configuration Objects %%%%%%

%% Carrier Configuration
% Configure the carriers for all of the gNBs with different physical layer cell identities.

cellIds = randperm(1008,numgNBs) - 1;

% Configure carrier properties
carrier = repmat(nrCarrierConfig,1,numgNBs);
for gNBIdx = 1:numgNBs
    carrier(gNBIdx).NCellID = cellIds(gNBIdx);
end
validateCarriers(carrier);

%% PRS Configuration
% Configure PRS parameters for all of the gNBs. Configure the PRS for different gNBs such that no overlap exists to avoid the problem of hearability.

% Slot offsets of different PRS signals
prsSlotOffsets = 0:2:(2*numgNBs - 1);
prsIDs = randperm(4096,numgNBs) - 1;    % AA TODO - DIFFERENZA IN 5G-NR TRA prsIDs e cellIds ??

% Configure PRS properties
prs = nrPRSConfig;
prs.PRSResourceSetPeriod = [10 0];
prs.PRSResourceOffset = 0;
prs.PRSResourceRepetition = 1;
prs.PRSResourceTimeGap = 1;
prs.MutingPattern1 = [];
prs.MutingPattern2 = [];
prs.NumRB = 52;     % defines the PRS spectral extension: each RB in 5G-NR consists of 12 subcarrriers ==> PRS spectral extension = 12*52 = 624 subcarriers
prs.RBOffset = 0;   % ?? 
prs.CombSize = 12;
prs.NumPRSSymbols = 12;
prs.SymbolStart = 0;
prs = repmat(prs,1,numgNBs);
for gNBIdx = 1:numgNBs
    prs(gNBIdx).PRSResourceOffset = prsSlotOffsets(gNBIdx);
    prs(gNBIdx).NPRSID = prsIDs(gNBIdx);
end

%% PDSCH Configuration
% Configure the physical downlink shared channel (PDSCH) parameters for all 
% of the gNBs for the data transmission. This example assumes that the data
% is transmitted from all of the gNBs that have a single transmission layer.
pdsch = nrPDSCHConfig;
pdsch.PRBSet = 0:51;
pdsch.SymbolAllocation = [0 14];
pdsch.DMRS.NumCDMGroupsWithoutData = 1;
pdsch = repmat(pdsch,1,numgNBs);
validateNumLayers(pdsch);


%% Path Loss Configuration - ( TODO EDIT/CHANGE FOR BETTER SCENARIO MODELING ) 
% Create a path loss configuration object for the 'UMa' (Urban Macrocell) scenario
plCfg = nrPathLossConfig;
plCfg.Scenario = 'Uma';
plCfg.BuildingHeight = 5; % [m] It is required for 'RMa' scenario
plCfg.StreetWidth = 5;  % [m] It is required for 'RMa' scenario
plCfg.EnvironmentHeight = 1; % [m] It is required for 'UMa' and 'UMi' scenario

% Specify the flag to configure the existence of the line of sight (LOS) path between each gNB and UE pair.
%los = [true false true true true]; % NB: 3 gNBs in los ==> 2 spatial unknowns retrievable ==> 2D positioning
los = ones(1,numgNBs);

if numel(los) ~= numgNBs
    error('nr5g:InvalidLOSLength',['Length of line of sight flag (' num2str(numel(los)) ...
        ') must be equal to the number of configured gNBs (' num2str(numgNBs) ').']);
end


%%%%%% Generate PRS and PDSCH Resources %%%%%%

%% Generate PRS and PDSCH resources corresponding to all of the gNBs 
% To improve the hearability of PRS signals, transmit PDSCH resources in the slots in which the PRS is not transmitted by any of the gNBs. 
% This example generates PRS and PDSCH resources with unit power (0 dB).

totSlots = nFrames*carrier(1).SlotsPerFrame;
prsGrid = cell(1,numgNBs);
dataGrid = cell(1,numgNBs);

for slotIdx = 0:totSlots-1
    [carrier(:).NSlot] = deal(slotIdx);
    [prsSym,prsInd] = deal(cell(1,numgNBs));
    for gNBIdx = 1:numgNBs
        % Create an empty resource grid spanning one slot in time domain
        slotGrid = nrResourceGrid(carrier(gNBIdx),1);

        % Generate PRS symbols and indices
        prsSym{gNBIdx} = nrPRS(carrier(gNBIdx),prs(gNBIdx));
        prsInd{gNBIdx} = nrPRSIndices(carrier(gNBIdx),prs(gNBIdx));

        % Map PRS resources to slot grid
        slotGrid(prsInd{gNBIdx}) = prsSym{gNBIdx};
        prsGrid{gNBIdx} = [prsGrid{gNBIdx} slotGrid];
    end
    % Transmit data in slots in which the PRS is not transmitted by any of
    % the gNBs (to control the hearability problem)
    for gNBIdx = 1:numgNBs
        dataSlotGrid = nrResourceGrid(carrier(gNBIdx),1);
        if all(cellfun(@isempty,prsInd))
            % Generate PDSCH indices
            [pdschInd,pdschInfo] = nrPDSCHIndices(carrier(gNBIdx),pdsch(gNBIdx));

            % Generate random data bits for transmission
            data = randi([0 1],pdschInfo.G,1);
            % Generate PDSCH symbols
            pdschSym = nrPDSCH(carrier(gNBIdx),pdsch(gNBIdx),data);

            % Generate demodulation reference signal (DM-RS) indices and symbols
            dmrsInd = nrPDSCHDMRSIndices(carrier(gNBIdx),pdsch(gNBIdx));
            dmrsSym = nrPDSCHDMRS(carrier(gNBIdx),pdsch(gNBIdx));

            % Map PDSCH and its associated DM-RS to slot grid
            dataSlotGrid(pdschInd) = pdschSym;
            dataSlotGrid(dmrsInd) = dmrsSym;
        end
        dataGrid{gNBIdx} = [dataGrid{gNBIdx} dataSlotGrid];
    end
end

% Plot carrier grid
plotGrid(prsGrid,dataGrid);


%%%%%% Perform OFDM Modulation %%%%%%

%% Perform orthogonal frequency-division multiplexing (OFDM) modulation of the signals at each gNB.

% Perform OFDM modulation of PRS and data signal at each gNB
txWaveform = cell(1,numgNBs);
for waveIdx = 1:numgNBs
    carrier(waveIdx).NSlot = 0;

    txWaveform{waveIdx} = nrOFDMModulate(carrier(waveIdx), prsGrid{waveIdx} + dataGrid{waveIdx});
end


% Compute OFDM information using first carrier, assuming all carriers are
% at same sampling rate
ofdmInfo = nrOFDMInfo(carrier(1));


% FIXME: improve plotting function
plotWaveform(txWaveform{1}, ofdmInfo.SampleRate, "Tx (OFDM) waveform ", carrier(1), 1, numgNBs)


%%%%%% Add Signal Delays and Apply Path Loss %%%%%%

% Calculate the (TRUE) time delays from each gNB to UE by using the known gNB and UE positions. 

% This calculation uses the (TRUE) distance between the UE and gNB (ie "radius" variable), and the speed of propagation (speed of light). 

% Calculate the sample delay by using the sampling rate, ofdmInfo.SampleRate, and store it in sampleDelayTrue. 
% This example only applies the integer sample delay by rounding the actual delays to their nearest integers. 
% The example uses these variables to model the environment between gNBs and the UE. The information about these variables is not provided/known to the UE.

speedOfLight = physconst('LightSpeed'); % [m/s]
radius = cell(1, numgNBs);              % [m]   TRUE distance(gNB,UE) 
delay = cell(1, numgNBs);               % [s]   TRUE TOAs       
sampleDelayTrue = zeros(1,numgNBs);     % [samples]

for  gNBIdx = 1:numgNBs
    radius{gNBIdx} = sqrt( (gNBPos{gNBIdx}(1) - UEPos(1))^2 + (gNBPos{gNBIdx}(2) - UEPos(2))^2); % [m]
    delay{gNBIdx} = radius{gNBIdx}/speedOfLight;                          % [s]
    sampleDelayTrue(gNBIdx) = round(delay{gNBIdx}*ofdmInfo.SampleRate);   % Delay in [samples]        ( ofdmInfo.SampleRate u.m: [samples/s] )
end

% Model the received signal at the UE by delaying each gNB transmission according to the values in sampleDelayTrue 
% and by attenuating the received signal from each gNB. 
% Ensure that all of the waveforms are of the same length by padding the received waveform from each gNB with relevant number of zeros.

rxWaveform = zeros(length(txWaveform{1}) + max(sampleDelayTrue),1); 
rx = cell(1,numgNBs);
% PathLoss = cell(1,numgNBs); %% TODO: store here the PL from each gNB for better debug  
for gNBIdx = 1: numgNBs

    % Calculate path loss for each gNB and UE pair
    PLdB = nrPathLoss(plCfg,fc,los(gNBIdx),[gNBPos{gNBIdx}(:);0],[UEPos(:);0]);
    if PLdB < 0 || isnan(PLdB) || isinf(PLdB)
        error('nr5g:invalidPL',"Computed path loss (" + num2str(PLdB) + ...
            ") is invalid. Try changing the UE or gNB positions, or path loss configuration.");
    end
    PL = 10^(PLdB/10);

    % Add delay, pad, and attenuate the waveforms ( by sqrt(PL) ) 
    rx{gNBIdx} = [zeros(sampleDelayTrue(gNBIdx),1); txWaveform{gNBIdx}; ...
                zeros(max(sampleDelayTrue)-sampleDelayTrue(gNBIdx),1)]/sqrt(PL);

    % Sum waveforms from all gNBs
    rxWaveform = rxWaveform + rx{gNBIdx};
end


% FIXME: improve plotting function
figure;
plot((1:length(rxWaveform))',real(rxWaveform));
title('real part of Rx waveform')



%%%%%% TOA Estimation (from 5G-NR PRS signals) %%%%%%

% cellsToBeDetected = min(3,numgNBs); % FIXME: ok?? (original matlab code issue ??)
% if cellsToBeDetected < 3 || cellsToBeDetected > numgNBs % meaning?? cellsToBeDetected > numgNBs condition cannot be valid ever!
%     error('nr5g:InvalidNumDetgNBs',['The number of detected gNBs (' num2str(x) ...
%         ') must be greater than or equal to 3 and less than or equal to the total number of gNBs (' num2str(numgNBs) ').']);
% end

if (length(los(los==1))>=3) % OK ?? not applying any selection of best (correlation) gNB and assuming los vector is known ?? 
    cellsToBeDetected = length(los(los==1));
else
     error('nr5g:InvalidNumDetgNBs',['The number of detected gNBs (' num2str(x) ...
         ') must be greater than or equal to 3 and less than or equal to the total number of gNBs (' num2str(numgNBs) ').']);
end

corr = cell(1,numgNBs);
sampleDelayEst = zeros(1, numgNBs); % ESTIMATED TOAs [samples]

maxCorr = zeros(1,numgNBs);
offsets = zeros(1,numgNBs);
nrCheckIsValid = zeros(1,numgNBs);

for gNBIdx = 1:numgNBs
    %[~,mag] = nrTimingEstimate(carrier(gNBIdx),rxWaveform,prsGrid{gNBIdx}); 
    [offset,mag] = nrTimingEstimate(carrier(gNBIdx),rxWaveform,prsGrid{gNBIdx}); % offset u.m. = [samples]
    offsets(gNBIdx) = offset;

    % Extract correlation data samples spanning about 1/14 ms (??) for normal
    % cyclic prefix and about 1/12 ms (??) for extended cyclic prefix (this
    % truncation is to ignore noisy side lobe peaks in correlation outcome)
    corr{gNBIdx} = mag(1:(ofdmInfo.Nfft*carrier(1).SubcarrierSpacing/15));

    % Delay estimate is at point of maximum correlation
    maxCorr(gNBIdx) = max(corr{gNBIdx});

    sampleDelayEst(gNBIdx) = find(corr{gNBIdx} == maxCorr(gNBIdx),1)-1;   % [samples]
    nrCheckIsValid(gNBIdx) = offset == sampleDelayEst(gNBIdx);
end

% Get detected gNB numbers based on correlation outcome
[~,detectedgNBs] = sort(maxCorr,'descend');
detectedgNBs = detectedgNBs(1:cellsToBeDetected);
 
% Plot PRS correlation results
plotPRSCorr(carrier,corr,ofdmInfo.SampleRate);



%%%%%% UE Position Estimation Using OTDOA Technique %%%%%%

% Calculate TDOA or RSTD values for all pairs of gNBs using the estimated TOA values 
% by using the getRSTDValues function. Each RSTD value corresponding to a pair of gNBs 
% can result from the UE being located at any position on a hyperbola with foci located at these gNBs.


% Compute RSTD values (ie TDOAs) for multilateration or trilateration 
if ~validationMode % validationMode=true
   rstdVals = getRSTDValues(sampleDelayEst,ofdmInfo.SampleRate); % [s] 
else % validationMode=true
   %% NB: zeros (TDOAs) matrix of size=numgNBsxnumgNBs ( in case of gNBs equidistant from UE ) 
   rstdVals = getRSTDValues(sampleDelayTrue,ofdmInfo.SampleRate); % [s] 
end

% Plot gNB and UE positions
txCellIDs = [carrier(:).NCellID];
cellIdx = 1;

curveX = {};
curveY = {};


for jj = detectedgNBs(1) % AA NB: Assuming FIRST detected gNB as reference/serving gNB
    for ii = detectedgNBs(2:end)

        rstd = rstdVals(ii,jj)*speedOfLight; % Delay distance [m] ( Range-Difference Measurement [m] )
       

        % Establish gNBs for which delay distance is applicable by
        % examining detected cell identities
        txi = find(txCellIDs == carrier(ii).NCellID);
        txj = find(txCellIDs == carrier(jj).NCellID);


        if (~isempty(txi) && ~isempty(txj))
            

            % Get x- and y- cartesian coordinates of hyperbola curve and store them
            %[x,y] = getRSTDCurve(gNBPos{txi},gNBPos{txj},rstd);
            [x,y, delta, phi, r, hk] = getRSTDCurve_improved(gNBPos{txi},gNBPos{txj},rstd);

            if isreal(x) && isreal(y)
                curveX{1,cellIdx} = x;
                curveY{1,cellIdx} = y;

                % Get the gNB numbers corresponding to the current hyperbola curve
                gNBNums{cellIdx} = [jj ii];
                cellIdx = cellIdx + 1;
            else
                warning("Due to the error in timing estimations, " + ...
                    "hyperbola curve cannot be deduced for the combination of gNB" + jj + " and gNB" + ii+". So, ignoring this gNB pair for the estimation of UE position.");
            end
        end
    end
end
numHyperbolaCurves = numel(curveX);
if numHyperbolaCurves < 2
    error('nr5g:InsufficientCurves',"The number of hyperbola curves ("+numHyperbolaCurves+") is not sufficient for the estimation of UE position in 2-D. Try with more number of gNBs.");
end




% Estimate UE position from hyperbola curves using
% |getEstimatedUEPosition|. This function computes the point of
% intersections of all hyperbola curves and does the average of those
% points to get the UE position.
estimatedUEPos = getEstimatedUEPosition(curveX,curveY);


% Compute positioning estimation error
EstimationErr = norm(UEPos-estimatedUEPos); % [m]

disp(['True UE Position [m]     : [' num2str(UEPos(1)) ' ' num2str(UEPos(2)) ']']);
disp(newline)


disp(['Estimated UE Position (OTDOA GEOMETRIC INTERSECTION method)      : [' num2str(estimatedUEPos(1)) ' ' num2str(estimatedUEPos(2)) ']' 13 ...
      'UE Position Estimation Error: ' num2str(EstimationErr) ' meters']);


% Plot UE, gNB positions, and hyperbola curves
gNBsToPlot = unique([gNBNums{:}],'stable');

plotPositionsAndHyperbolaCurves(gNBPos,UEPos,gNBsToPlot,curveX,curveY,gNBNums,estimatedUEPos);




diary off
movefile(logFileName, outputDirectory);
    














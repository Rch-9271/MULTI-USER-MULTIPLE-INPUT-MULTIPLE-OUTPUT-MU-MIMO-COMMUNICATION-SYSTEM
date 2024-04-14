clc;                              % Clear worksapace
clear all;
close all;
s = rng(67);                        % Set RNG state for repeatability

%% DEFINE SYSTEM PARAMETERS FOR MULTI-USER MULTIPLE-INPUT MULTIPLE-OUTPUT (MU-MIMO) COMMUNICATION SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining Multi-user system with single/multiple streams per user
config.numUsers = 4;                   % Number of users  in the system.
config.numSTSVec = [3 2 1 2];          % Vector indicating the Number of independent data streams per user
config.numSTS = sum(config.numSTSVec);    % Total number of spatial streams across all users. It's the sum of prm.numSTSVec. Must be a power of 2
config.numTx = config.numSTS * 8;         % Number of BS transmit antennas (power of 2)
config.numRx = config.numSTSVec * 4;      % Number of receive antennas, per user.It's a vector with each element indicating the number of antennas for each user. (any >= numSTSVec)

% Each user has the same modulation
config.bitsPerSubCarrier = 4;          % Number of bits per subcarrier used for modulation  (2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM)
config.numDataSymbols = 10;            % Number of OFDM data symbols

% MS positions: assumes BS at origin
maxRange = 1000; %   Angles specified as [azimuth;elevation] degrees
                 %   az in range [-180 180], el in range [-90 90], e.g. [45;0]
                 % all MSs within 1000 meters of BS
config.mobileRanges = randi([1 maxRange], 1, config.numUsers);
config.mobileAngles = [rand(1, config.numUsers) * 360 - 180; ...
                    rand(1, config.numUsers) * 180 - 90];

config.fc = 28e9;                      % 28 GHz system
config.chanSRate = 100e6;              % Channel sampling rate, 100 Msps
config.ChanType = 'Scattering';        % Channel options: 'Scattering', 'MIMO'
config.NFig = 8;                       % Noise figure (increase to worsen, 5-10 dB)
config.nRays = 500;                    % Number of rays for Frf, Fbb partitioning

%% DEFINE OFDM MODULATION PARAMETERS USEDS FOR THE SYSTEM   %%%%%%%%%%%%%%%

config.FFTLength = 256;
config.CyclicPrefixLength = 64;
config.numCarriers = 234;              % Number of carries
config.NullCarrierIndices = [1:7 129 256-5:256]'; % Guards and DC
config.PilotCarrierIndices = [26 54 90 118 140 168 204 232]';
nonDataIdx = [config.NullCarrierIndices; config.PilotCarrierIndices];
config.CarriersLocations = setdiff((1:config.FFTLength)', sort(nonDataIdx));

numSTS = config.numSTS;
numTx = config.numTx;
numRx = config.numRx;
numSTSVec = config.numSTSVec;
codeRate = 1/3;                     % same code rate per user
numTails = 6;                       % number of termination tail bits
config.numFrmBits = numSTSVec.* (config.numDataSymbols * config.numCarriers * ...
                 config.bitsPerSubCarrier * codeRate) - numTails;
config.modMode = 2 ^ config.bitsPerSubCarrier; % Modulation order
% Account for channel filter delay
numPadSym = 3;                      % number of symbols to zeropad
config.numPadZeros = numPadSym* (config.FFTLength + config.CyclicPrefixLength);

%% DEFINE TRANSMITE AND RECEIVE ARRAYS AND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% POSITIONAL PARAMETERS FOR THE SYSTEM %%%%%%%%%%%%%%%%%%%%

config.cLight = physconst('LightSpeed');
config.lambda = config.cLight / config.fc;

% Get transmit and receive array information
[isTxURA, expFactorTx, isRxURA, expFactorRx] = helperArrayInfo(config, true);

% Transmit antenna array definition
%   Array locations and angles
config.posTx = [0; 0; 0];              % BS/Transmit array position, [x; y; z], meters
if isTxURA
    % Uniform Rectangular array
    txArray = phased.PartitionedArray(...
        'Array', phased.URA([expFactorTx numSTS], 0.5 * config.lambda), ...
        'SubarraySelection', ones(numSTS, numTx), 'SubarraySteering', 'Custom');
else
    % Uniform Linear array
    txArray = phased.ULA(numTx, 'ElementSpacing', 0.5 * config.lambda, ...
        'Element', phased.IsotropicAntennaElement('BackBaffled', false));
end
config.posTxElem = getElementPosition(txArray) / config.lambda;

spLoss = zeros(config.numUsers, 1);
config.posRx = zeros(3, config.numUsers);
for uIdx = 1:config.numUsers

    % Receive arrays
    if isRxURA(uIdx)
        % Uniform Rectangular array
        rxarray = phased.PartitionedArray(...
            'Array', phased.URA([expFactorRx(uIdx) numSTSVec(uIdx)], ...
            0.5*config.lambda), 'SubarraySelection', ones(numSTSVec(uIdx), ...
            numRx(uIdx)), 'SubarraySteering', 'Custom');
        config.posRxElem = getElementPosition(rxarray) / config.lambda;
    else
        if numRx(uIdx) > 1
            % Uniform Linear array
            rxarray = phased.ULA(numRx(uIdx), ...
                'ElementSpacing', 0.5 * config.lambda, ...
                'Element', phased.IsotropicAntennaElement);
            config.posRxElem = getElementPosition(rxarray)/config.lambda;
        else
            rxarray = phased.IsotropicAntennaElement;
            config.posRxElem = [0; 0; 0]; % LCS
        end
    end

    % Mobile positions
    [xRx, yRx, zRx] = sph2cart(deg2rad(config.mobileAngles(1, uIdx)), ...
                             deg2rad(config.mobileAngles(2, uIdx)), ...
                             config.mobileRanges(uIdx));
    config.posRx(:,uIdx) = [xRx; yRx; zRx];
    [toRxRange, toRxAng] = rangeangle(config.posTx, config.posRx(:, uIdx));
    spLoss(uIdx) = fspl(toRxRange, config.lambda);
end


%% CHANNEL STATE INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the preamble signal
config.numSTS = numTx;                 % set to numTx to sound out all channels
preambleSig = helperGenPreamble(config);

% Transmit preamble over channel
config.numSTS = numSTS; % keep same array config for channel
[rxPreSig, chanDelay] = helperApplyMUChannel(preambleSig, config, spLoss);

% Channel state information feedback
hDp = cell(config.numUsers, 1);
config.numSTS = numTx;                 % set to numTx to estimate all links
for uIdx = 1:config.numUsers

    % Front-end amplifier gain and thermal noise
    rxPreAmp = phased.ReceiverPreamp( ...
        'Gain', spLoss(uIdx), ...    % account for path loss
        'NoiseFigure', config.NFig,'ReferenceTemperature',290, ...
        'SampleRate', config.chanSRate);
    rxPreSigAmp = rxPreAmp(rxPreSig{uIdx});
    %   scale power for used sub-carriers
    rxPreSigAmp = rxPreSigAmp * (sqrt(config.FFTLength - ...
        length(config.NullCarrierIndices)) / config.FFTLength);

    % OFDM demodulation
    rxOFDM = ofdmdemod(rxPreSigAmp(chanDelay(uIdx) + 1: ...
        end - (config.numPadZeros-chanDelay(uIdx)), :), config.FFTLength, ...
        config.CyclicPrefixLength, config.CyclicPrefixLength, ...
        config.NullCarrierIndices, config.PilotCarrierIndices);

    % Channel estimation from preamble
    %       numCarr, numTx, numRx
    hDp{uIdx} = helperMIMOChannelEstimate(rxOFDM(:, 1:numTx, :), config);

end


%% HYBRID BEAMFORMING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the hybrid weights on the transmit side
if config.numUsers == 1
    % Single-user OMP
    %   Spread rays in [az;el]=[-180:180;-90:90] 3D space, equal spacing
    %   txang = [-180:360/prm.nRays:180; -90:180/prm.nRays:90];
    txang = [rand(1, config.nRays) * 360 - 180; rand(1, config.nRays) * 180 - 90]; % random
    At = steervec(config.posTxElem, txang);

    Fbb = complex(zeros(config.numCarriers, numSTS, numSTS));
    Frf = complex(zeros(config.numCarriers, numSTS, numTx));
    for carrIdx = 1:config.numCarriers
        [Fbb(carrIdx, :, :), Frf(carrIdx, :, :)] = helperOMPTransmitWeights( ...
            permute(hDp{1}(carrIdx, :, :), [2 3 1]), numSTS, numSTS, At);
    end
    v = Fbb;    % set the baseband precoder (Fbb)
    % Frf is same across subcarriers for flat channels
    mFrf = permute(mean(Frf, 1), [2 3 1]);

else
    % Multi-user Joint Spatial Division Multiplexing
    [Fbb, mFrf] = helperJSDMTransmitWeights(hDp, config);

    % Multi-user baseband precoding
    %   Pack the per user CSI into a matrix (block diagonal)
    steeringMatrix = zeros(config.numCarriers, sum(numSTSVec), sum(numSTSVec));
    for uIdx = 1:config.numUsers
        stsIdx = sum(numSTSVec(1:uIdx-1)) + (1:numSTSVec(uIdx));
        steeringMatrix(:, stsIdx, stsIdx) = Fbb{uIdx};  % Nst-by-Nsts-by-Nsts
    end
    v = permute(steeringMatrix, [1 3 2]);

end

% Transmit array pattern plots
if isTxURA
    % URA element response for the first subcarrier
    pattern(txArray, config.fc, -180:180, -90:90, 'Type', 'efield', ...
            'ElementWeights', mFrf.' * squeeze(v(1, :, :)), ...
            'PropagationSpeed', config.cLight);
else % ULA
    % Array response for first subcarrier
    wts = mFrf.' * squeeze(v(1, :, :));
    pattern(txArray, config.fc, -180:180, -90:90, 'Type', 'efield', ...
            'Weights', wts(:, 1), 'PropagationSpeed', config.cLight);
end

config.numSTS = numSTS;                 % revert back for data transmission


%% Data Transmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convolutional encoder
encoder = comm.ConvolutionalEncoder( ...
    'TrellisStructure',poly2trellis(7,[133 171 165]), ...
    'TerminationMethod','Terminated');

% Bits to QAM symbol mapping
modRQAM = comm.RectangularQAMModulator( ...
    'ModulationOrder',config.modMode,'BitInput',true, ...
    'NormalizationMethod','Average power');

txDataBits = cell(config.numUsers, 1);
gridData = complex(zeros(config.numCarriers,config.numDataSymbols,numSTS));
for uIdx = 1:config.numUsers
    % Generate mapped symbols from bits per user
    txDataBits{uIdx} = randi([0,1],config.numFrmBits(uIdx),1);
    encodedBits = encoder(txDataBits{uIdx});
    mappedSym = modRQAM(encodedBits);

    % Map to layers: per user, per symbol, per data stream
    stsIdx = sum(numSTSVec(1:(uIdx-1)))+(1:numSTSVec(uIdx));
    gridData(:,:,stsIdx) = reshape(mappedSym,config.numCarriers, ...
        config.numDataSymbols,numSTSVec(uIdx));
end

% Apply precoding weights to the subcarriers, assuming perfect feedback
preData = complex(zeros(config.numCarriers,config.numDataSymbols,numSTS));
for symIdx = 1:config.numDataSymbols
    for carrIdx = 1:config.numCarriers
        Q = squeeze(v(carrIdx,:,:));
        normQ = Q * sqrt(numTx)/norm(Q,'fro');
        preData(carrIdx,symIdx,:) = squeeze(gridData(carrIdx,symIdx,:)).' ...
            * normQ;
    end
end

% Multi-antenna pilots
pilots = helperGenPilots(config.numDataSymbols,numSTS);

% OFDM modulation of the data
txOFDM = ofdmmod(preData,config.FFTLength,config.CyclicPrefixLength,...
                 config.NullCarrierIndices,config.PilotCarrierIndices,pilots);
%   scale power for used sub-carriers
txOFDM = txOFDM * (config.FFTLength/ ...
    sqrt((config.FFTLength-length(config.NullCarrierIndices))));

% Generate preamble with the feedback weights and prepend to data
preambleSigD = helperGenPreamble(config,v);
txSigSTS = [preambleSigD;txOFDM];

% RF beamforming: Apply Frf to the digital signal
%   Each antenna element is connected to each data stream
txSig = txSigSTS*mFrf;


%% Signal Propagation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply a spatially defined channel to the transmit signal
[rxSig,chanDelay] = helperApplyMUChannel(txSig,config,spLoss,preambleSig);


%% Receive Amplification and Signal Recovery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hfig = figure('Name','Equalized symbol constellation per stream');
scFact = ((config.FFTLength-length(config.NullCarrierIndices))...
         /config.FFTLength^2)/numTx;
nVar = noisepow(config.chanSRate,config.NFig,290)/scFact;
demodRQAM = comm.RectangularQAMDemodulator( ...
    'ModulationOrder',config.modMode,'BitOutput',true, ...
    'DecisionMethod','Approximate log-likelihood ratio', ...
    'NormalizationMethod','Average power','Variance',nVar);
decoder = comm.ViterbiDecoder('InputFormat','Unquantized', ...
    'TrellisStructure',poly2trellis(7, [133 171 165]), ...
    'TerminationMethod','Terminated','OutputDataType','double');

%% Show simple summary

fprintf('\n----- Execution Summary -----\n');

fprintf('Number of users: %d\n', config.numUsers);
disp(['Number of independent data streams per user: [' num2str(config.numSTSVec) ']']);
fprintf('numSTS: %d\n', config.numSTS);
fprintf('Number of BS transmit antennas: %d\n', config.numTx);
disp(['Number of receive antennas, per user: [' num2str(config.numRx) ']']);
fprintf('Number of bits per sub carrier: %d\n', config.bitsPerSubCarrier);

if config.bitsPerSubCarrier == 2
    fprintf('Modulation: QPSK\n')
elseif config.bitsPerSubCarrier == 4
    fprintf('Modulation: 16QAM\n')
elseif config.bitsPerSubCarrier == 6
    fprintf('Modulation: 64QAM\n')
elseif config.bitsPerSubCarrier == 8
    fprintf('Modulation: 256QAM\n')
end

fprintf('Number of OFDM data symbols: %d\n', config.numDataSymbols);
fprintf('maxRange: %d\n', maxRange);
fprintf('Frequency: %d\n', config.fc)
fprintf('Maximum Sample Rate: %d\n', config.chanSRate)
fprintf('Channel type: %s\n', config.ChanType)
fprintf('Noise figure: %d\n', config.NFig)
fprintf('Number of rays for Frf: %d\n\n', config.nRays)

fprintf('\n-----------------------------\n');

%% Show RMS EVM and BER for each user %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for uIdx = 1:config.numUsers
    stsU = numSTSVec(uIdx);
    stsIdx = sum(numSTSVec(1:(uIdx-1)))+(1:stsU);

    % Front-end amplifier gain and thermal noise
    rxPreAmp = phased.ReceiverPreamp( ...
        'Gain', spLoss(uIdx), ...        % account for path loss
        'NoiseFigure', config.NFig,'ReferenceTemperature',290, ...
        'SampleRate', config.chanSRate);
    rxSigAmp = rxPreAmp(rxSig{uIdx});

    % Scale power for occupied sub-carriers
    rxSigAmp = rxSigAmp*(sqrt(config.FFTLength-length(config.NullCarrierIndices)) ...
        /config.FFTLength);

    % OFDM demodulation
    rxOFDM = ofdmdemod(rxSigAmp(chanDelay(uIdx)+1: ...
        end-(config.numPadZeros-chanDelay(uIdx)),:),config.FFTLength, ...
        config.CyclicPrefixLength,config.CyclicPrefixLength, ...
        config.NullCarrierIndices,config.PilotCarrierIndices);

    % Channel estimation from the mapped preamble
    hD = helperMIMOChannelEstimate(rxOFDM(:,1:numSTS,:),config);

    % MIMO equalization
    %   Index into streams for the user of interest
    [rxEq,CSI] = helperMIMOEqualize(rxOFDM(:,numSTS+1:end,:),hD(:,stsIdx,:));

    % Soft demodulation
    rxSymbs = rxEq(:)/sqrt(numTx);
    rxLLRBits = demodRQAM(rxSymbs);

    % Apply CSI prior to decoding
    rxLLRtmp = reshape(rxLLRBits,config.bitsPerSubCarrier,[], ...
        config.numDataSymbols,stsU);
    csitmp = reshape(CSI,1,[],1,numSTSVec(uIdx));
    rxScaledLLR = rxLLRtmp.*csitmp;

    % Soft-input channel decoding
    rxDecoded = decoder(rxScaledLLR(:));

    % Decoded received bits
    rxBits = rxDecoded(1:config.numFrmBits(uIdx));

    % Plot equalized symbols for all streams per user
    scaler = ceil(max(abs([real(rxSymbs(:)); imag(rxSymbs(:))])));
    for i = 1:stsU
        subplot(config.numUsers, max(numSTSVec), (uIdx-1)*max(numSTSVec)+i);
        plot(reshape(rxEq(:,:,i)/sqrt(numTx), [], 1), '.');
        axis square
        xlim(gca,[-scaler scaler]);
        ylim(gca,[-scaler scaler]);
        title(['U ' num2str(uIdx) ', DS ' num2str(i)]);
        grid on;
    end

    % Compute and display the EVM
    evm = comm.EVM('Normalization','Average constellation power', ...
        'ReferenceSignalSource','Estimated from reference constellation', ...
        'ReferenceConstellation',constellation(demodRQAM));
    rmsEVM = evm(rxSymbs);
    disp(['User ' num2str(uIdx)]);
    disp(['  RMS EVM (%) = ' num2str(rmsEVM)]);

    % Compute and display bit error rate
    ber = comm.ErrorRate;
    measures = ber(txDataBits{uIdx},rxBits);
    fprintf('  BER = %.5f; No. of Bits = %d; No. of errors = %d\n', ...
        measures(1),measures(3),measures(2));
end

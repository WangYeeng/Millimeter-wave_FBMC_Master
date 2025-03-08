%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2024 by Ying Wang.   
%==========================================================================
clear;clc;close all;
load('Results\BER_FBMC_Alamouti_1.mat');
BER_FBMC_Alamouti_1 = BER_FBMC_Alamouti;
clear BER_FBMC_Alamouti;
load('Results\BER_FBMC_AWGN.mat');
BER_FBMC_AWGN       = BER_FBMC_Alamouti;
clear BER_FBMC_Alamouti;
rng(0);
%% Parameters
NrRepetitions       = 100;               % Number of Monte Carlo repetitions
M_SNR_dB            = -10:5:15;          % Signal-to-noise ratio
QAM_ModulationOrder = 4;                 % Modulation order, 4,16,64,...
L                   = 256;               % Number of subcarriers
K                   = 16;                % Number of FBMC symbols per Block = Spreading Length
Kall                = K*3+3;             % Total Number of FBMC symbols (3 guard symbols!)

SubcarrierSpacing = 120e3;               % For 'CDL-D_3ns' with L = 256, the Opt.Subcarrier Spacing is 120 kHZ. 
CarrierFrequency  = 60e9;
SamplingRate      = SubcarrierSpacing*L*15;  % Sampling rate (must be larger than the subcarrier spacing times L)

Dopplershift      = 1;
%% FBMC Object
FBMC = Modulation.FBMC(...
    L,...                                 % Number of subcarriers
    Kall,...                              % Number of FBMC symbols
    SubcarrierSpacing,...                 % Subcarrier spacing (Hz)
    SamplingRate,...                      % Sampling rate (Samples/s)
    CarrierFrequency,...                  % Intermediate frequency first subcarrier (Hz)
    false,...                             % Transmit real valued signal
    'PHYDYAS-OQAM',...                    % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    4, ...                                % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                % Initial phase shift
    true ...                              % Polyphase implementation
    );

%% Modulation Object
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
% For ML Detection
ML_MapIndex1 = reshape(repmat((1:QAM.ModulationOrder),QAM.ModulationOrder,1),1,QAM.ModulationOrder^2);
ML_MapIndex2 = reshape(repmat((1:QAM.ModulationOrder).',1,QAM.ModulationOrder),1,QAM.ModulationOrder^2);
ML_Mapping = QAM.SymbolMapping([ML_MapIndex1;ML_MapIndex2]);
%% Channel
if Dopplershift
    Velocity_kmh      = 3;
    PowerDelayProfile = 'CDL-D_8ns';                               % Power delay profile

    ChannelModel = Channel.FastFading(...
                SamplingRate,...                                   % Sampling rate (Samples/s)
                PowerDelayProfile,...                              % Power delay profile, either string or vector
                FBMC.Nr.SamplesTotal,...                           % Number of total samples
                Velocity_kmh/3.6*CarrierFrequency/2.998e8,...      % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
                'Jakes',...                                        % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
                3,...                                              % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum
                2,...                                              % Number of transmit antennas
                2,...                                              % Number of receive antennas
                false ...                                          % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
                );   
end

%% Alamouti Object
Alamouti = Modulation.SpaceCoding(...
    'Alamouti2x1',...                     % Space time coding method
    1 ...                                 % Frequency spreading = 0; time spreading = 1
    );

% Get Coding matrix
FBMC.SetNrMCSymbols(K);
C = PrecodingMatrix.GetPrecodingMatrix(FBMC,1);
CBlock = [zeros(L,L*K/2)  ;...  
                 C         ]; 
FBMC.SetNrMCSymbols(Kall);            
% Power Normalization so that Transmitted Power is the same in OFDM and FBMC
FBMCPowerNormalization = sqrt(2*(K+1)/K); % 2: FBMC-OQAM structure and (K+1)/K because of zero guard symbol

%% Pilot Matrix: 0=DataSymbol, 1=PilotSymbol, -1=ZeroSymbol;
PilotMatrixBlockAntenna1 = zeros(L,K/2);
PilotMatrixBlockAntenna1(2:6:end,1:8:end)=1; 
PilotMatrixBlockAntenna1(5:6:end,5:8:end)=1;
PilotMatrixBlockAntenna1(2:6:end,2:8:end)=-1; 
PilotMatrixBlockAntenna1(5:6:end,6:8:end)=-1;
PilotMatrixBlockAntenna2 = PilotMatrixBlockAntenna1*(-1);
PilotMatrixBlockAntenna(:,:,1) = PilotMatrixBlockAntenna1;
PilotMatrixBlockAntenna(:,:,2) = PilotMatrixBlockAntenna2;

NrPilots = sum(PilotMatrixBlockAntenna1(:)==1);
NrDataSymbols = sum(PilotMatrixBlockAntenna1(:)==0);
NrTransmittedSymbols = length(PilotMatrixBlockAntenna1(:));
%%
BER_FBMC_Alamouti       = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_SM_ZeroForcing = nan(length(M_SNR_dB),NrRepetitions);
BER_FBMC_SM_ML          = nan(length(M_SNR_dB),NrRepetitions);

%% Simulate Over Different Channel Realizations
tic;
for i_rep = 1:NrRepetitions  
%% Generate Data and Pilots
% Pilot Symbols: The pilot symbol power is by a factor of two higher because (pilots for the other antenna) are zero!
x_PilotAntenna1 = QAM.SymbolMapping(randi(QAM.ModulationOrder,NrPilots,1));
x_PilotAntenna2 = QAM.SymbolMapping(randi(QAM.ModulationOrder,NrPilots,1));
x_PilotAntenna1 = x_PilotAntenna1./abs(x_PilotAntenna1)*sqrt(2);
x_PilotAntenna2 = x_PilotAntenna2./abs(x_PilotAntenna2)*sqrt(2);

% Binary Data Stream
BinaryDataStream_Alamouti   = randi([0 1],NrDataSymbols*log2(QAM.ModulationOrder),1);
BinaryDataStream_SMAntenna1 = randi([0 1],NrDataSymbols*log2(QAM.ModulationOrder),1); %Spatial Multiplexing
BinaryDataStream_SMAntenna2 = randi([0 1],NrDataSymbols*log2(QAM.ModulationOrder),1); %Spatial Multiplexing

% Transmitted Alamouti Symbols
x_Alamouti = nan(L,K/2);
x_Alamouti(PilotMatrixBlockAntenna1==0) = QAM.Bit2Symbol(BinaryDataStream_Alamouti);
x_Alamouti_Coded = Alamouti.Encoder(x_Alamouti);
x_Alamouti_Coded(PilotMatrixBlockAntenna==1)=[x_PilotAntenna1;x_PilotAntenna2];
x_Alamouti_Coded(PilotMatrixBlockAntenna==-1)=0;
x_Alamouti_Coded_Antenna1 = x_Alamouti_Coded(:,:,1);
x_Alamouti_Coded_Antenna2 = x_Alamouti_Coded(:,:,2);

% Transmitted Spatial Multiplexed Symbols
x_SM_Antenna1 = nan(L,K/2);
x_SM_Antenna1(PilotMatrixBlockAntenna1==0)  = QAM.Bit2Symbol(BinaryDataStream_SMAntenna1);
x_SM_Antenna1(PilotMatrixBlockAntenna1==1)  = x_PilotAntenna1;
x_SM_Antenna1(PilotMatrixBlockAntenna1==-1) = 0;
x_SM_Antenna2 = nan(L,K/2);
x_SM_Antenna2(PilotMatrixBlockAntenna2==0)  = QAM.Bit2Symbol(BinaryDataStream_SMAntenna2);
x_SM_Antenna2(PilotMatrixBlockAntenna2==1)  = x_PilotAntenna2;
x_SM_Antenna2(PilotMatrixBlockAntenna2==-1) = 0;

% Data symbols of the the first and second block (chosen randomly to keep it simple)
x_Block1_Antenna1 = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));
x_Block1_Antenna2 = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));
x_Block2_Antenna1 = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));
x_Block2_Antenna2 = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));

% Account for pilots in block 1 and 3
x_Block1_Antenna1(PilotMatrixBlockAntenna1==1)  = x_Block1_Antenna1(PilotMatrixBlockAntenna1==1)./abs(x_Block1_Antenna1(PilotMatrixBlockAntenna1==1))*sqrt(2);
x_Block1_Antenna2(PilotMatrixBlockAntenna2==1)  = x_Block1_Antenna2(PilotMatrixBlockAntenna2==1)./abs(x_Block1_Antenna2(PilotMatrixBlockAntenna2==1))*sqrt(2);
x_Block1_Antenna1(PilotMatrixBlockAntenna1==-1) = 0;
x_Block1_Antenna2(PilotMatrixBlockAntenna2==-1) = 0;
x_Block2_Antenna1(PilotMatrixBlockAntenna1==1)  = x_Block2_Antenna1(PilotMatrixBlockAntenna1==1)./abs(x_Block2_Antenna1(PilotMatrixBlockAntenna1==1))*sqrt(2);
x_Block2_Antenna2(PilotMatrixBlockAntenna2==1)  = x_Block2_Antenna2(PilotMatrixBlockAntenna2==1)./abs(x_Block2_Antenna2(PilotMatrixBlockAntenna2==1))*sqrt(2);
x_Block2_Antenna1(PilotMatrixBlockAntenna1==-1) = 0;
x_Block2_Antenna2(PilotMatrixBlockAntenna2==-1) = 0;


%% Transmitted Signal
TransmittedSymbols_FBMC_Alamouti_Antenna1 = reshape(CBlock*[x_Block1_Antenna1(:) x_Alamouti_Coded_Antenna1(:) x_Block2_Antenna1(:)],L,K*3+3);
TransmittedSymbols_FBMC_Alamouti_Antenna2 = reshape(CBlock*[x_Block1_Antenna2(:) x_Alamouti_Coded_Antenna2(:) x_Block2_Antenna2(:)],L,K*3+3);
s_FBMC_Alamouti_Antenna1 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_Alamouti_Antenna1);
s_FBMC_Alamouti_Antenna2 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_Alamouti_Antenna2);

TransmittedSymbols_FBMC_SM_Antenna1 = reshape(CBlock*[x_Block1_Antenna1(:) x_SM_Antenna1(:) x_Block2_Antenna1(:)],L,K*3+3);
TransmittedSymbols_FBMC_SM_Antenna2 = reshape(CBlock*[x_Block1_Antenna2(:) x_SM_Antenna2(:) x_Block2_Antenna2(:)],L,K*3+3);
s_FBMC_SM_Antenna1 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_SM_Antenna1);
s_FBMC_SM_Antenna2 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_SM_Antenna2);

if Dopplershift
    %% Fast fading channel model
    ChannelModel.NewRealization;
    h11 = ChannelModel.GetConvolutionMatrix{1,1};
    h12 = ChannelModel.GetConvolutionMatrix{1,2};
    h21 = ChannelModel.GetConvolutionMatrix{2,1};
    h22 = ChannelModel.GetConvolutionMatrix{2,2};

else
    %% Lens Anetnna Array
    % Saleh-Valenzuela channel model
    M1 = 1;                 % number of vertical anetnnas
    M2 = 2;                 % number of horizental antennas
    M_BS = M1*M2;
    Q  = 2;                 % number of single antenna users
    Np = 10; ArrayGeometry = 'UPA';
    H_sv=zeros(M1*M2,Q); 
    for q=1:Q
       H_sv(:,q) = Channel.Saleh_Valenzuelachannel(M1,M2,Np,ArrayGeometry);
    end
    % the beamspace channel
    H = Channel.UPA_LensArray(M1,M2,Q,H_sv);
    % H = Channel.ULA_LensArray(M_BS,Q,1,H_sv);
    h11  = H(1,1);
    h12  = H(1,2);
    h21  = H(2,1);
    h22  = H(2,2);
end

%% Simulate Over Different Noise Values
for i_SNR = 1:length(M_SNR_dB)
SNR_dB = M_SNR_dB(i_SNR); 

Pn = SamplingRate/(SubcarrierSpacing*L)*10^(-SNR_dB/10); 
noise_Antenna1 = sqrt(Pn/2)*(randn(size(s_FBMC_Alamouti_Antenna1))+1j*randn(size(s_FBMC_Alamouti_Antenna1)));
noise_Antenna2 = sqrt(Pn/2)*(randn(size(s_FBMC_Alamouti_Antenna1))+1j*randn(size(s_FBMC_Alamouti_Antenna1)));

%% Received Signal
r_FBMC_Alamouti_Antenna1 = h11*s_FBMC_Alamouti_Antenna1+h12*s_FBMC_Alamouti_Antenna2+noise_Antenna1;

r_FBMC_SM_Antenna1 = h11*s_FBMC_SM_Antenna1+h12*s_FBMC_SM_Antenna2+noise_Antenna1;
r_FBMC_SM_Antenna2 = h21*s_FBMC_SM_Antenna1+h22*s_FBMC_SM_Antenna2+noise_Antenna2;

%% Demodulation
y_FBMC_Alamouti_3Blocks = FBMC.Demodulation(r_FBMC_Alamouti_Antenna1);
y_FBMC_Alamouti = reshape(CBlock'*reshape(y_FBMC_Alamouti_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);

y_FBMC_SM_Antenna1_3Blocks = FBMC.Demodulation(r_FBMC_SM_Antenna1);
y_FBMC_SM_Antenna1 = reshape(CBlock'*reshape(y_FBMC_SM_Antenna1_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);
y_FBMC_SM_Antenna2_3Blocks = FBMC.Demodulation(r_FBMC_SM_Antenna2);
y_FBMC_SM_Antenna2 = reshape(CBlock'*reshape(y_FBMC_SM_Antenna2_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);

%% Channel Estimation
% Note that the noise power is also increasd=> The same SNR (ignoring zero guard)
EstimatedChannel_FBMC_Alamouti(:,:,1) = mean(y_FBMC_Alamouti(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1)*ones(L,K/2);
EstimatedChannel_FBMC_Alamouti(:,:,2) = mean(y_FBMC_Alamouti(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2)*ones(L,K/2);

H11_Est_FBMC = mean(y_FBMC_SM_Antenna1(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
H21_Est_FBMC = mean(y_FBMC_SM_Antenna2(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
H12_Est_FBMC = mean(y_FBMC_SM_Antenna1(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);
H22_Est_FBMC = mean(y_FBMC_SM_Antenna2(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);

H_Est_FBMC =[H11_Est_FBMC,H12_Est_FBMC;H21_Est_FBMC,H22_Est_FBMC];
%% Data Detection
%Estimated Symbols
x_est_FBMC_Alamouti       = Alamouti.Decoder(y_FBMC_Alamouti,EstimatedChannel_FBMC_Alamouti*sqrt(2));
x_est_FBMC_SM_ZeroForcing = reshape((pinv(H_Est_FBMC)*[y_FBMC_SM_Antenna1(:) y_FBMC_SM_Antenna2(:)].').',L,K/2,2);
% ML Detection
y_FBMC_SM_Temp = repmat(reshape([y_FBMC_SM_Antenna1(:).';y_FBMC_SM_Antenna2(:).'],2,1,[]),1,QAM.ModulationOrder^2,1);
[~,indexMin] = min(sum(abs(y_FBMC_SM_Temp-repmat(H_Est_FBMC*ML_Mapping,1,1,L*K/2)).^2,1),[],2);
x_est_FBMC_SM_ML = reshape(ML_Mapping(:,indexMin(:)).',L,K/2,2);

% Symbols To Bit
DetectedBitStream_FBMC_Alamouti       = QAM.Symbol2Bit(x_est_FBMC_Alamouti(PilotMatrixBlockAntenna1==0));
DetectedBitStream_FBMC_SM_ZeroForcing = QAM.Symbol2Bit(x_est_FBMC_SM_ZeroForcing(PilotMatrixBlockAntenna==0));
DetectedBitStream_FBMC_SM_ML          = QAM.Symbol2Bit(x_est_FBMC_SM_ML(PilotMatrixBlockAntenna==0));
% Bit Error Ratio
BER_FBMC_Alamouti(i_SNR,i_rep)       = mean(BinaryDataStream_Alamouti~=DetectedBitStream_FBMC_Alamouti);
BER_FBMC_SM_ZeroForcing(i_SNR,i_rep) = mean([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]~=DetectedBitStream_FBMC_SM_ZeroForcing);
BER_FBMC_SM_ML(i_SNR,i_rep)   = mean([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]~=DetectedBitStream_FBMC_SM_ML);

end

TimePassed = toc;
if mod(i_rep,10)==0
disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes']);
end
end
%% Plot Results
MarkerSize = 10;
Linewidth  = 1.2;
figure(1);
semilogy(M_SNR_dB,mean(BER_FBMC_SM_ZeroForcing,2),'-o','Color',[0 0 1]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
hold on;grid on;
semilogy(M_SNR_dB,mean(BER_FBMC_SM_ML,2),         '-d','Color',[0 1 0]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
% semilogy(M_SNR_dB,mean(BER_FBMC_Alamouti_1,2),    '-^','Color',[1 1 1]*0.45,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
semilogy(M_SNR_dB,mean(BER_FBMC_Alamouti,2),      '-^','Color',[1 0 0]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
if Dopplershift
    semilogy(M_SNR_dB,mean(BER_FBMC_AWGN,2),            '-.','Color',[1 1 1]*0.25,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
    legend('MIMO-ZF','MIMO-ML','MISO-Alamouti2x1 code','Theory, AWGN','Location','southwest');
    ylim([1e-3 1e0]);
    title('CDL-D, 3km/h@60GHz, $N_P = 10$','Interpreter','latex');
else
    semilogy(M_SNR_dB,mean(BER_FBMC_AWGN,2),            '-.','Color',[1 1 1]*0.25,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
    legend('MIMO-ZF','MIMO-ML','MISO-Alamouti2x1 code 1 Path','MISO-Alamouti2x1 code 10 Paths','Theory, AWGN','Location','southwest');
    ylim([1e-3 1e0]);
    title('Lens Antenna Array, UPA, $N_P = 10$','Interpreter','latex');
end
xlabel('Signal-to-Noise Ratio, $P_s/P_n$ (dB)','Interpreter','latex');
ylabel('Bit Error Ratio');
set(gca,'FontName','Times New Roman','FontSize',12,'GridLineWidth',1.2,'LooseInset', [0,0,0,0]);
%% Plot Additional Information
fprintf('=============================\n');
fprintf('========= Data Rate =========\n');
fprintf('FBMC       |%7.2f Mbit/s  | \n', length([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]) / (FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols)/1e6   );
fprintf('=============================\n');

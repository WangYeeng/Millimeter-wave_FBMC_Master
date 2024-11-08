%==========================================================================
% Ying Wang, wangyingstu@163.com
% (c) 2024 by Ying Wang.   
%==========================================================================
clear;clc;close all;
%% Parameters
NrRepetitions       = 100;               % Number of Monte Carlo repetitions
M_SNR_dB            = -5:4:20;           % Signal-to-noise ratio
QAM_ModulationOrder = 4;                 % Modulation order, 4,16,64,...
L                   = 128;               % Number of subcarriers
K                   = 32;                % Number of FBMC symbols per Block = Spreading Length
Kall                = K*3+3;             % Total Number of FBMC symbols (3 guard symbols!)

SubcarrierSpacing = 90e3;                 
CarrierFrequency  = 60e9;
SamplingRate      = SubcarrierSpacing*L*15;  % Sampling rate (must be larger than the subcarrier spacing times L)

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
NMSE_FBMC          = nan(length(M_SNR_dB),NrRepetitions);

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
TransmittedSymbols_FBMC_SM_Antenna1 = reshape(CBlock*[x_Block1_Antenna1(:) x_SM_Antenna1(:) x_Block2_Antenna1(:)],L,K*3+3);
TransmittedSymbols_FBMC_SM_Antenna2 = reshape(CBlock*[x_Block1_Antenna2(:) x_SM_Antenna2(:) x_Block2_Antenna2(:)],L,K*3+3);
s_FBMC_SM_Antenna1 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_SM_Antenna1);
s_FBMC_SM_Antenna2 = FBMCPowerNormalization/sqrt(2)*FBMC.Modulation(TransmittedSymbols_FBMC_SM_Antenna2);

%% 透镜阵列矩阵
% Saleh-Valenzuela channel model
M1 = 1; % number of vertical anetnnas
M2 = 2; % number of horizental antennas
Q  = 2;  % number of single antenna users
Np = 10; ArrayGeometry = 'UPA';
H_sv=zeros(M1*M2,Q); % the beamspace channel 
for q=1:Q
  H_sv(:,q) = Channel.Saleh_Valenzuelachannel(M1,M2,Np,ArrayGeometry);
end
H = Channel.UPA_lensarray(M1,M2,Q,H_sv);
h11  = H(1,1);
h12  = H(1,2);
h21  = H(2,1);
h22  = H(2,2);
%% Simulate Over Different Noise Values
for i_SNR = 1:length(M_SNR_dB)
SNR_dB = M_SNR_dB(i_SNR); 

Pn = SamplingRate/(SubcarrierSpacing*L)*10^(-SNR_dB/10); 
noise_Antenna1 = sqrt(Pn/2)*(randn(size(s_FBMC_SM_Antenna1))+1j*randn(size(s_FBMC_SM_Antenna1)));
noise_Antenna2 = sqrt(Pn/2)*(randn(size(s_FBMC_SM_Antenna2))+1j*randn(size(s_FBMC_SM_Antenna2)));

%% Received Signal
r_FBMC_SM_Antenna1 = h11*s_FBMC_SM_Antenna1+h12*s_FBMC_SM_Antenna2+noise_Antenna1;
r_FBMC_SM_Antenna2 = h21*s_FBMC_SM_Antenna1+h22*s_FBMC_SM_Antenna2+noise_Antenna2;

%% Demodulation
y_FBMC_SM_Antenna1_3Blocks = FBMC.Demodulation(r_FBMC_SM_Antenna1);
y_FBMC_SM_Antenna1 = reshape(CBlock'*reshape(y_FBMC_SM_Antenna1_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);
y_FBMC_SM_Antenna2_3Blocks = FBMC.Demodulation(r_FBMC_SM_Antenna2);
y_FBMC_SM_Antenna2 = reshape(CBlock'*reshape(y_FBMC_SM_Antenna2_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);

%% Channel Estimation
H11_Est_FBMC = mean(y_FBMC_SM_Antenna1(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
H21_Est_FBMC = mean(y_FBMC_SM_Antenna2(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
H12_Est_FBMC = mean(y_FBMC_SM_Antenna1(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);
H22_Est_FBMC = mean(y_FBMC_SM_Antenna2(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);

H_Est_FBMC = [H11_Est_FBMC,H12_Est_FBMC;H21_Est_FBMC,H22_Est_FBMC];
H_Perfect  = [h11,h12;h21,h22];

NMSE_FBMC(i_SNR,i_rep) = 10*log10((norm(abs(H_Est_FBMC)-abs(H_Perfect)))^2/(norm(H_Perfect))^2);

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
semilogy(M_SNR_dB,mean(NMSE_FBMC,2),'-o','Color',[0 0 1]*0.75,'MarkerSize',MarkerSize,'Linewidth',Linewidth);
hold on;grid on;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('NMSE (dB)');
legend('UPA','Location','southwest');
set(gca,'FontName','Times New Roman','FontSize',12,'GridLineWidth',1.2,'LooseInset', [0,0,0,0]);


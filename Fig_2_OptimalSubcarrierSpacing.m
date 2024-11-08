% =====================================================================    
% (c) 2024 by Ying Wang, wangyingstu@163.com
% =====================================================================    
% This script calculates the Signal-to-Interference Ratio (SIR) for
% different subcarrier spacings.
% =====================================================================
clear; close all; clc;rng(0);
%% Parameters
L                      = 32;
K_FBMC                 = 7;
M_SubcarrierSpacing    = (10:2:150)*1e3;
CarrierFrequency       = 60e9;
PrototypeFilter        = 'PHYDYAS';           % Prototype filter.
PowerDelayProfile      = 'ExtendedVehicularA';
Velocity_kmh           = 40;  
% PowerDelayProfile      = 'CDL-D_5ns';
% Velocity_kmh           = 3;
%% variables space
F_Theory_opt_PHYDYAS    = nan(1,length(M_SubcarrierSpacing));
M_RealSubcarrierSpacing = nan(1,length(M_SubcarrierSpacing));
SIR_Theory_Flat         = nan(1,length(M_SubcarrierSpacing));
SIR_FBMC                = nan(1,length(M_SubcarrierSpacing));
%% Start Simulation
for i_subcarrier = 1:length(M_SubcarrierSpacing)
    SubcarrierSpacing = M_SubcarrierSpacing(i_subcarrier);
    %% 设置采样率
    if  SubcarrierSpacing <= 25e3
        Samplingrate = SubcarrierSpacing*L*5;
    elseif (25e3<SubcarrierSpacing)&&SubcarrierSpacing <= 80e3
        Samplingrate = SubcarrierSpacing*L*15;
    else
        Samplingrate = SubcarrierSpacing*L*12;
    end
   
    Multicarrier = Modulation.FBMC(...
                L,...                                           % Number subcarriers.
                K_FBMC,...                                      % Number FBMC symbols.
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                Samplingrate,...                                % Sampling rate (Samples/s)
                CarrierFrequency,...                            % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmit real valued signal
                [PrototypeFilter '-OQAM'],...                   % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
                4, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
                0, ...                                          % Initial phase shift
                true ...                                        % Polyphase implementation
                );
   TXMatrix = Multicarrier.GetTXMatrix;
   M_RealSubcarrierSpacing(i_subcarrier) = Multicarrier.PHY.SubcarrierSpacing;
    %% Channel Object
ChannelModel = Channel.FastFading(...
        Samplingrate,...                                        % Sampling rate (Samples/s)
        PowerDelayProfile,...                                   % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
        Multicarrier.Nr.SamplesTotal,...                        % Number of total samples
        Velocity_kmh/3.6*CarrierFrequency/2.998e8,...           % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
        'Jakes',...                                             % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
        200, ...                                                % Number multipath delays for WSS process
        1,...                                                   % Number of transmit antennas
        1,...                                                   % Number of receive antennas
        0 ...                                                   % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );
    R_vecH = ChannelModel.GetCorrelationMatrix;
    TimeOffset = round(ChannelModel.GetMeanDelay/ChannelModel.PHY.dt);
    RMSDelaySpread = ChannelModel.GetRmsDelaySpread;
    %% Channel property approximation for F
    F_Theory_opt_PHYDYAS(i_subcarrier)            = 0.91*sqrt(ChannelModel.PHY.MaximumDopplerShift/(sqrt(2)*RMSDelaySpread));
    %% Hypergeometric function approximation for flat fading
    P_Signal_Theory_Flat = hypergeom(1/2,[3/2 2], -(pi*ChannelModel.PHY.MaximumDopplerShift...
                                          /(Multicarrier.PHY.SubcarrierSpacing)).^2);
    P_Interference_Theory_Flat = 1-P_Signal_Theory_Flat;
    SIR_Theory_Flat(i_subcarrier)   = 10*log10(P_Signal_Theory_Flat/P_Interference_Theory_Flat);
    %% Simulation measured F
    [PSignal,PInterference] = Multicarrier.GetSignalAndInterferencePowerOQAM(...
                R_vecH,...
                eye(Multicarrier.Nr.Subcarriers*Multicarrier.Nr.MCSymbols),...
                TimeOffset,...
                ceil(Multicarrier.Nr.Subcarriers/2),...
                ceil(Multicarrier.Nr.MCSymbols/2));
    SIR_FBMC(i_subcarrier) = 10*log10(PSignal./PInterference);
    
    disp(i_subcarrier/length(M_SubcarrierSpacing));
end
F_Theory_opt_PHYDYAS  = mean(F_Theory_opt_PHYDYAS);
[SIR_max,Index_F_opt] = max(SIR_FBMC);
F_opt                 = M_SubcarrierSpacing(Index_F_opt);
%% Plot Results
Linewidth = 1.5;
plot(M_RealSubcarrierSpacing/1e3,SIR_Theory_Flat                ,'Color',[1 0 0]*0.75,'Linewidth',Linewidth);
hold on;grid on;
plot(M_RealSubcarrierSpacing/1e3,SIR_FBMC                       ,'Color',[0 0 1]*0.75,'Linewidth',Linewidth);
plot([F_Theory_opt_PHYDYAS/1e3 F_Theory_opt_PHYDYAS/1e3],[10 32],'-.','Color',[1 1 1]*0.25,'Linewidth',Linewidth);
plot([F_opt/1e3 F_opt/1e3],[10 32]                              ,'Color',[0 1 0]*0.65,'Linewidth',Linewidth);
xlabel('Subcarrier spacing (kHz)');
ylabel('Signal-to-Interference Ratio (dB)');
%set(gca,'Xtick',[0 50 60 68 100 150],'FontName','Times New Roman','FontSize',12,'GridLineWidth',1.2,'GridLineStyle','-.','LooseInset', [0,0,0,0]);
set(gca,'FontName','Times New Roman','FontSize',12,'GridLineWidth',1.2,'GridLineStyle','-.','LooseInset', [0,0,0,0]);


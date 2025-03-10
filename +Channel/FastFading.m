classdef FastFading < handle
    % =====================================================================       
    % The convolution is performed by the method
    % ".Convolution(.)". A new channel realization is obtained by 
    % ".NewRealization;".
    % Additionally, there exist some other useful methods, such as
    % ".GetRmsDelaySpread" or ".PlotPowerDelayProfile"
    % =====================================================================    
    
    
    properties (SetAccess = private)
        PHY
        Nr
        Implementation
        ImpulseResponse
    end
    
    methods
        function obj = FastFading(...
                SamplingRate,...                                   % Sampling rate (Samples/s)
                PowerDelayProfile,...                              % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
                SamplesTotal,...                                   % Number of total samples
                MaximumDopplerShift,...                            % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
                DopplerModel,...                                   % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
                Paths,...                                          % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum
                nTxAntennas,...                                    % Number of transmit antennas
                nRxAntennas,...                                    % Number of receive antennas
                WarningIfSampleRateDoesNotMatch ...                % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
                )
            
            % Initialize parameters
            obj.PHY.SamplingRate                    = SamplingRate;
            obj.Nr.SamplesTotal                     = SamplesTotal;
            obj.Nr.txAntennas                       = nTxAntennas;
            obj.Nr.rxAntennas                       = nRxAntennas;
            obj.PHY.MaximumDopplerShift             = MaximumDopplerShift;
            obj.Implementation.PowerDelayProfile    = PowerDelayProfile;
            
            obj.PHY.dt = 1/obj.PHY.SamplingRate;
            
            if ischar(PowerDelayProfile)
                % See for Pedestrian A ... ftp://www.3gpp.org/tsg_ran/TSG_RAN/TSGR_16/Docs/PDF/RP-020376.pdf
                if strcmp(PowerDelayProfile(1:3),'CDL')
                    Pos1 = strfind(PowerDelayProfile,'_');
                    Pos2 = strfind(PowerDelayProfile,'ns');
                    DesiredRMSdelaySpread = str2double(PowerDelayProfile(Pos1+1:Pos2-1))*10^-9;
                    PowerDelayProfile = PowerDelayProfile(1:5);
                end
                switch PowerDelayProfile
                    case 'Flat'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0;...                                           % Relative power in dB
                            0];                                             % Relative delay in ns
                    case 'AWGN'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0;...                                           % Relative power in dB
                            0];                                             % Relative delay in ns                        
                        if (obj.PHY.MaximumDopplerShift>0) && strcmp(PowerDelayProfile,'AWGN')
                            warning('The velocity is set to zero (AWGN channel!)');
                            obj.PHY.MaximumDopplerShift = 0;
                        end
                    case 'PedestrianA'
                        PowerDelayProfileINdB_DelayIndex = [...
                            0,-9.7,-19.2,-22.8;...                          % Relative power in dB
                            0,110e-9,190e-9,410e-9];                        % Relative delay in ns
                    case 'ExtendedVehicularA'
                      PowerDelayProfileINdB_DelayIndex = [...                      
                          0,-1.5,-1.4,-3.6,-0.6,-9.1,-7,-12,-16.9; 
                          0,30e-9,150e-9,310e-9,370e-9,710e-9,1090e-9,1730e-9,2510e-9];
                    case 'CDL-D'
                        PowerDelayProfileINdB_DelayIndex(1,:) = ...
                            [-0.2,-13.5,-18.8,-21,-22.8,-17.9,-20.1,-21.9,-22.9,-27.8,-23.6,-24.8,-30,-27.7];
                        PowerDelayProfileINdB_DelayIndex(2,:) = ...
                            DesiredRMSdelaySpread*[0.000,0.000,0.035,0.612,1.363,1.405,1.804,2.596,1.775,4.042,7.937,9.424,9.708,12.525];
                    case 'CDL-E'
                        PowerDelayProfileINdB_DelayIndex(1,:) = ...
                            [-0.03,-22.03,-15.8,-18.1,-19.8,-22.9,-22.4,-18.6,-20.8,-22.6,-22.3,-25.6,-20.2,-29.8,-29.2];
                        PowerDelayProfileINdB_DelayIndex(2,:) = ...
                            DesiredRMSdelaySpread*[0.0000,0.0000,0.5133,0.5440,0.5630,0.5440,0.7112,1.9092,1.9293,1.9589,2.6426,3.7136,5.4524,12.0034,20.6419];
                    otherwise
                        error('Power delay profile model not supported!');
                end
                IndexDelays = round(PowerDelayProfileINdB_DelayIndex(2,:)./obj.PHY.dt)+1;
                if  WarningIfSampleRateDoesNotMatch && (sum(abs(rem(PowerDelayProfileINdB_DelayIndex(2,:),obj.PHY.dt)))>0)
                    disp('Sampling rate does not match the predefined delays of the channel model!');
                    disp('Delay taps are changed according to (top = desired value; bottom = chosen value):');
                    disp([PowerDelayProfileINdB_DelayIndex(2,:);(IndexDelays-1)*obj.PHY.dt]);
                end
                PowerDelayProfileTemp = zeros(length(IndexDelays),max(IndexDelays));
                for i_delay = 1:length(IndexDelays)
                    PowerDelayProfileTemp(i_delay,IndexDelays(i_delay)) = 10.^(PowerDelayProfileINdB_DelayIndex(1,i_delay)/10);
                end
                obj.PHY.PowerDelayProfile = sum(PowerDelayProfileTemp,1);
                obj.PHY.DesiredPowerDelayProfiledB = PowerDelayProfileINdB_DelayIndex; 
            else
                obj.PHY.PowerDelayProfile = PowerDelayProfile;
                
                obj.PHY.DesiredPowerDelayProfiledB(1,:) = 10*log10(PowerDelayProfile);
                obj.PHY.DesiredPowerDelayProfiledB(2,:) = (0:size(obj.PHY.DesiredPowerDelayProfiledB,2)-1)*obj.PHY.dt;                
            end
            obj.Implementation.PowerDelayProfileNormalized = obj.PHY.PowerDelayProfile.'/sum(obj.PHY.PowerDelayProfile);
            
            obj.Implementation.IndexDelayTaps = find(obj.PHY.PowerDelayProfile');
            
            NrOnes = sum(obj.Nr.SamplesTotal:-1:(obj.Nr.SamplesTotal-length(obj.PHY.PowerDelayProfile)+1));
            obj.Implementation.MappingConvolutionMatrix = nan(NrOnes,2);
            for i_NrTaps = 1:length(obj.PHY.PowerDelayProfile)
                obj.Implementation.MappingConvolutionMatrix(sum(obj.Nr.SamplesTotal:-1:obj.Nr.SamplesTotal-i_NrTaps+2)+(1:obj.Nr.SamplesTotal-i_NrTaps+1),:) = [-i_NrTaps+1+(i_NrTaps:obj.Nr.SamplesTotal).' (i_NrTaps:obj.Nr.SamplesTotal).'];
            end
           
            % Faster mapping because zeros are not mapped
            obj.Implementation.MappingConvolutionMatrixFast = obj.Implementation.MappingConvolutionMatrix;
            for i_NrTaps = find(obj.PHY.PowerDelayProfile==0)
                obj.Implementation.MappingConvolutionMatrixFast(sum(obj.Nr.SamplesTotal:-1:obj.Nr.SamplesTotal-i_NrTaps+2)+(1:obj.Nr.SamplesTotal-i_NrTaps+1),:) = -1;
            end  
            obj.Implementation.MappingConvolutionMatrixFast(obj.Implementation.MappingConvolutionMatrixFast(:,1)==-1,:)=[];
            
            if obj.PHY.MaximumDopplerShift/(obj.PHY.SamplingRate/obj.Nr.SamplesTotal)<=0.5 && strcmp(DopplerModel(1:min(8,end)),'Discrete') && obj.PHY.MaximumDopplerShift>0
                warning('Discrete Doppler spectrum: The velocity is so low, that it is set to zero. Change to a "Jakes" or "Uniform" spectrum or increase the number of symbols in time!');
                obj.PHY.MaximumDopplerShift = 0;
            end            

            if obj.PHY.MaximumDopplerShift>0
                obj.PHY.DopplerModel        = DopplerModel;
                if strcmp(obj.PHY.DopplerModel(1:min(8,end)),'Discrete')
                    obj.Implementation.UseDiscreteDopplerSpectrum = true;
                    % Calculate frequency spacing between discrete Doppler shifts
                    df = obj.PHY.SamplingRate/obj.Nr.SamplesTotal;

                    NrDopplerShifts = ceil(obj.PHY.MaximumDopplerShift/df);

                    IntervalPoints = df*((-NrDopplerShifts-1:NrDopplerShifts)+0.5);
                    IntervalPoints(IntervalPoints<=-obj.PHY.MaximumDopplerShift)=-obj.PHY.MaximumDopplerShift;
                    IntervalPoints(IntervalPoints>=obj.PHY.MaximumDopplerShift)=obj.PHY.MaximumDopplerShift;

                    switch obj.PHY.DopplerModel
                        case 'Discrete-Jakes'
                            DiscreteDopplerSpectrum = asin(IntervalPoints(2:end)/obj.PHY.MaximumDopplerShift)-asin(IntervalPoints(1:end-1)/obj.PHY.MaximumDopplerShift); % calculate the integral of the continous Doppler spectrum over a certain interval
                        case 'Discrete-Uniform'
                            DiscreteDopplerSpectrum = IntervalPoints(2:end)-IntervalPoints(1:end-1); % calculate the integral of the continous Doppler spectrum over a certain interval
                        otherwise
                            error('Doppler spectrum not supported');
                    end    

                    % Normalize
                    DiscreteDopplerSpectrum = DiscreteDopplerSpectrum./sum(DiscreteDopplerSpectrum);

                    obj.Implementation.DiscreteDopplerSpectrum = repmat(DiscreteDopplerSpectrum.',1,length(obj.Implementation.IndexDelayTaps));            

                else
                    obj.Nr.Paths                = Paths;                    
                    obj.Implementation.UseDiscreteDopplerSpectrum = false;
                end
            end
  
            obj.NewRealization;
            
            obj.Implementation.CancelElementsConvolutionMatrix = ones(obj.Nr.SamplesTotal,length(obj.PHY.PowerDelayProfile))==1;
            for i_NrTaps = 1:length(obj.PHY.PowerDelayProfile)-1
                obj.Implementation.CancelElementsConvolutionMatrix(i_NrTaps,(end-length(obj.PHY.PowerDelayProfile)+1+i_NrTaps):end) = false;
            end
            obj.Implementation.CancelElementsConvolutionMatrixFast = obj.Implementation.CancelElementsConvolutionMatrix(:,obj.Implementation.IndexDelayTaps);
                 
        end
        
        function NewRealization( obj )
            % obtain new realization of the channel

            if strcmp(obj.Implementation.PowerDelayProfile, 'AWGN')
                obj.ImpulseResponse(1,1,:,:) = ones( obj.Nr.rxAntennas, obj.Nr.txAntennas );
            else
                if obj.PHY.MaximumDopplerShift>0
                    obj.ImpulseResponse = zeros( obj.Nr.SamplesTotal, size(obj.Implementation.PowerDelayProfileNormalized,1), obj.Nr.txAntennas, obj.Nr.txAntennas  );
                    ImpulseResponseTemp = zeros(length(obj.PHY.PowerDelayProfile),obj.Nr.SamplesTotal);
                    if obj.Implementation.UseDiscreteDopplerSpectrum % Faster calulation using IFFT. Underlying assumption: discrete Doppler shifts => valid only for a high number of samples/Doppler shift
                        for iTx = 1:obj.Nr.txAntennas
                            for iRx = 1:obj.Nr.rxAntennas
                                NrDopplerShifts = (size(obj.Implementation.DiscreteDopplerSpectrum,1)-1)/2;

                                GaussUncorr1 = obj.Nr.SamplesTotal/sqrt(2)*(randn(NrDopplerShifts+1,length(obj.Implementation.IndexDelayTaps))+1j*randn(NrDopplerShifts+1,length(obj.Implementation.IndexDelayTaps)));
                                GaussUncorr2 = obj.Nr.SamplesTotal/sqrt(2)*(randn(NrDopplerShifts,length(obj.Implementation.IndexDelayTaps))+1j*randn(NrDopplerShifts,length(obj.Implementation.IndexDelayTaps)));    
                     
                                GaussUncorr1 = GaussUncorr1.*repmat(sqrt(obj.Implementation.PowerDelayProfileNormalized(obj.Implementation.IndexDelayTaps).'),size(GaussUncorr1,1),1);
                                GaussUncorr2 = GaussUncorr2.*repmat(sqrt(obj.Implementation.PowerDelayProfileNormalized(obj.Implementation.IndexDelayTaps).'),size(GaussUncorr2,1),1);
                                
                                obj.ImpulseResponse(:,obj.Implementation.IndexDelayTaps,iRx,iTx) = ifft(...
                                     [sqrt(obj.Implementation.DiscreteDopplerSpectrum(NrDopplerShifts+1:end,:)).*GaussUncorr1;...
                                      zeros(obj.Nr.SamplesTotal-NrDopplerShifts*2-1,length(obj.Implementation.IndexDelayTaps));...
                                      sqrt(obj.Implementation.DiscreteDopplerSpectrum(1:NrDopplerShifts,:)).*GaussUncorr2]);    
                                                                
                  
                            end
                        end
                    else % This is the true Doppler spectum (continious and not discrete!)
                        for iTx = 1:obj.Nr.txAntennas
                            for iRx = 1:obj.Nr.rxAntennas
                                switch obj.PHY.DopplerModel
                                    case 'Jakes'
                                        DopplerShifts = cos(rand([length(obj.Implementation.IndexDelayTaps) 1 obj.Nr.Paths])*2*pi)*obj.PHY.MaximumDopplerShift;
                                    case 'Uniform'
                                        DopplerShifts = 2*(rand([length(obj.Implementation.IndexDelayTaps) 1 obj.Nr.Paths])-0.5)*obj.PHY.MaximumDopplerShift;
                                    otherwise
                                        error('Doppler spectrum not supported');
                                end
                                RandomPhase = rand([length(obj.Implementation.IndexDelayTaps) 1 obj.Nr.Paths]);
                                t = 0:obj.PHY.dt:obj.PHY.dt*(obj.Nr.SamplesTotal-1);
                                ImpulseResponseTemp(obj.Implementation.IndexDelayTaps,:) = sum(exp(1j*2*pi*(bsxfun(@plus,RandomPhase,bsxfun(@times,DopplerShifts,t)))),3)./sqrt(obj.Nr.Paths);

                                obj.ImpulseResponse(:,:,iRx,iTx) = bsxfun(@times,sqrt(obj.Implementation.PowerDelayProfileNormalized), ImpulseResponseTemp).';
                            end
                        end
                    end
                else
                    obj.ImpulseResponse = zeros(1,size(obj.Implementation.PowerDelayProfileNormalized,1), obj.Nr.txAntennas, obj.Nr.txAntennas  );
                    for iTx = 1:obj.Nr.txAntennas
                        for iRx = 1:obj.Nr.rxAntennas
                            obj.ImpulseResponse(1,:,iRx,iTx) = (1/sqrt(2)*sqrt(obj.Implementation.PowerDelayProfileNormalized).*(randn(size(obj.Implementation.PowerDelayProfileNormalized))+1j*randn(size(obj.Implementation.PowerDelayProfileNormalized)))).';
                        end
                    end
                end
            end
        end 
        
        
        function convolvedSignal = Convolution(obj, signal)
            % Convolution of the signal with the time variant impulse response
            
            N = size(signal,1);
            if obj.PHY.MaximumDopplerShift>0
                convolvedSignal = zeros( N, obj.Nr.rxAntennas);
                ConvolutionMatrix   = obj.GetConvolutionMatrix;
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        convolvedSignal(:,iRx) = convolvedSignal(:,iRx) + ConvolutionMatrix{iRx,iTx}(1:N,1:N)*signal(:,iTx);
                    end
                end
            else
                convolvedSignal = zeros( N + size(obj.Implementation.PowerDelayProfileNormalized,1)-1, obj.Nr.rxAntennas);
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        convolvedSignal(:,iRx) = convolvedSignal(:,iRx) + conv( signal(:,iTx), obj.ImpulseResponse(1,:,iRx,iTx) );
                    end
                end
                convolvedSignal = convolvedSignal(1:N,:);% remove tails
            end
        end
        
        function ConvolutionMatrix = GetConvolutionMatrix(obj)
            % returns the time-variant convolution matrix
            
            ConvolutionMatrix = cell( obj.Nr.rxAntennas, obj.Nr.txAntennas );
            if obj.PHY.MaximumDopplerShift>0
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        ImpulseResponseTemp = obj.ImpulseResponse(:,obj.Implementation.IndexDelayTaps,iRx,iTx);
                        ConvolutionMatrix{iRx,iTx} = sparse(obj.Implementation.MappingConvolutionMatrixFast(:,2),obj.Implementation.MappingConvolutionMatrixFast(:,1),ImpulseResponseTemp(obj.Implementation.CancelElementsConvolutionMatrixFast),obj.Nr.SamplesTotal,obj.Nr.SamplesTotal);
                    end
                end
            else
                for iTx = 1:obj.Nr.txAntennas
                    for iRx = 1:obj.Nr.rxAntennas
                        ImpulseResponseTemp = obj.ImpulseResponse(ones(obj.Nr.SamplesTotal,1),obj.Implementation.IndexDelayTaps,iRx,iTx);
                        ConvolutionMatrix{iRx,iTx} = sparse(obj.Implementation.MappingConvolutionMatrixFast(:,2),obj.Implementation.MappingConvolutionMatrixFast(:,1),ImpulseResponseTemp(obj.Implementation.CancelElementsConvolutionMatrixFast),obj.Nr.SamplesTotal,obj.Nr.SamplesTotal);
                    end
                end
            end
        end
        
        function ChannelTransferFunction = GetTransferFunction(obj, TimePos, FFTSize, ActiveSubcarrier)
            % Calculates the channel transfer function at certain time
            % positions. The first argument is a vector and represent the time
            % positions. The second argument represents the FFT size of the
            % underlying modulation format. The third argument (if
            % specified) determines which subcarriers are active (only 
            % these frequency positions are returned).


            ChannelTransferFunction = zeros( FFTSize, length(TimePos), obj.Nr.rxAntennas, obj.Nr.txAntennas );
            if obj.PHY.MaximumDopplerShift==0
                TimePos = ones(length(TimePos),1);
            end
            for iTx = 1:obj.Nr.txAntennas
                for iRx = 1:obj.Nr.rxAntennas
                    ImpulseTemp = [obj.ImpulseResponse(TimePos,:,iRx,iTx).';zeros(FFTSize-size(obj.ImpulseResponse,2),length(TimePos))];
                    ChannelTransferFunction(:,:,iRx,iTx) = fft(ImpulseTemp,[],1);
                end
            end  
            if exist('ActiveSubcarrier','var')
                ChannelTransferFunction = ChannelTransferFunction(ActiveSubcarrier,:,:,:);
            end
        end     
        
        function [TimeCorrelation,Time]=GetTimeCorrelation(obj)
            % returns the time-autocorrelation function of the channel.

            if obj.PHY.MaximumDopplerShift>0
                if obj.Implementation.UseDiscreteDopplerSpectrum
                    % Can be implemented similar as the frequency
                    % correlation. Needs to be done when time
                    warning('Exact solution only for a continuous Doppler spectrum.');
                end            
                Time = -obj.PHY.dt*(obj.Nr.SamplesTotal-1):obj.PHY.dt:obj.PHY.dt*(obj.Nr.SamplesTotal-1);
                switch obj.PHY.DopplerModel
                    case 'Jakes'
                        TimeCorrelation = besselj(0,pi*2*obj.PHY.MaximumDopplerShift*Time);
                    case 'Uniform'
                        TimeCorrelation = sinc(2*obj.PHY.MaximumDopplerShift*Time);
                end
            else
                TimeCorrelation = ones(1,obj.Nr.SamplesTotal*2-1);
            end
        end
        
        function [FrequencyCorrelation,Frequency] = GetFrequencyCorrelation(obj)
            % returns the frequency-autocorrelation function of the channel.
            
            df = 1/(obj.Nr.SamplesTotal*obj.PHY.dt);
            
            FrequencyCorrelation = fft([obj.Implementation.PowerDelayProfileNormalized;zeros(obj.Nr.SamplesTotal-length(obj.Implementation.PowerDelayProfileNormalized),1)]);
            FrequencyCorrelation = circshift(FrequencyCorrelation,[ceil(obj.Nr.SamplesTotal/2) 1]);
            
            Frequency =((1:obj.Nr.SamplesTotal)-ceil(obj.Nr.SamplesTotal/2)-1)*df;
        end
        
        function MeanDelay = GetMeanDelay(obj)
            % returns the mean delay of the channel.
            
            Tau = (0:(length(obj.Implementation.PowerDelayProfileNormalized))-1)*obj.PHY.dt;
            MeanDelay = sum(Tau.*obj.Implementation.PowerDelayProfileNormalized.');
        end
        function RmsDelaySpread = GetRmsDelaySpread(obj)
            % returns the root mean squared doppler spread of the channel
            
            Tau = (0:(length(obj.Implementation.PowerDelayProfileNormalized))-1)*obj.PHY.dt;
            MeanDelay = obj.GetMeanDelay;
            RmsDelaySpread = sqrt(sum(Tau.^2.*obj.Implementation.PowerDelayProfileNormalized.')-MeanDelay.^2);
        end
        function CorrelationMatrix = GetCorrelationMatrix(obj)
            % returns the correlation matrix of the vectorized convolution 
            % matrix. Let r=H*s be the channel description with H beeing
            % a time-variant convolution matrix. Then, this method
            % returns R_vecH=E{H(:)H(:)'}.

            IndexTimeCorrelation = bsxfun(@plus,obj.Nr.SamplesTotal+(0:(obj.Nr.SamplesTotal-1)).',(0:-1:-(obj.Nr.SamplesTotal-1)));
            TimeCorrelation =  obj.GetTimeCorrelation;
            CorrelationMatrixOrdered = (...
                repmat(TimeCorrelation(IndexTimeCorrelation),[size(obj.Implementation.PowerDelayProfileNormalized,1) 1]).*...
                kron(sparse(obj.Implementation.PowerDelayProfileNormalized),ones(size(IndexTimeCorrelation))));
            IndexCorrMatrixConv = bsxfun(@plus,(1:(obj.Nr.SamplesTotal+1):obj.Nr.SamplesTotal^2).',(0:1:(size(obj.ImpulseResponse,2)-1)));
            %           CorrelationMatrixOld = sparse(...
            %             repmat(reshape(IndexCorrMatrixConv,[],1),obj.Nr.SamplesTotal,1),...
            %             reshape(kron(IndexCorrMatrixConv.',ones(size(IndexCorrMatrixConv,1),1)),[],1),...
            %             CorrelationMatrixOrdered(:));
            % Memory efficient solution
            SplittingFactor = ceil(numel(IndexCorrMatrixConv)*obj.Nr.SamplesTotal*8/1e9);
            while mod(obj.Nr.SamplesTotal/SplittingFactor,1)
                SplittingFactor = SplittingFactor+1;
            end
            if SplittingFactor>1
                CorrelationMatrix = sparse(obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1,obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1);
                Row = repmat(reshape(IndexCorrMatrixConv,[],1),obj.Nr.SamplesTotal/SplittingFactor,1);
                for i_split = 1:SplittingFactor
                    Index1Temp = (1:obj.Nr.SamplesTotal/SplittingFactor)+(i_split-1)*obj.Nr.SamplesTotal/SplittingFactor;
                    Index2Temp = (1:numel(CorrelationMatrixOrdered)/SplittingFactor)+(i_split-1)*numel(CorrelationMatrixOrdered)/SplittingFactor;
                    
                    Column = reshape(kron(IndexCorrMatrixConv(Index1Temp,:).',ones(size(IndexCorrMatrixConv,1),1)),[],1);
                    CorrelationMatrix = CorrelationMatrix+...
                        sparse(Row,Column,CorrelationMatrixOrdered(Index2Temp),...
                        obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1,...
                        obj.Nr.SamplesTotal^2+size(IndexCorrMatrixConv,2)-1);
                end
            else
                CorrelationMatrix = sparse(...
                    repmat(reshape(IndexCorrMatrixConv,[],1),obj.Nr.SamplesTotal,1),...
                    reshape(kron(IndexCorrMatrixConv.',ones(size(IndexCorrMatrixConv,1),1)),[],1),...
                    CorrelationMatrixOrdered(:));
            end
            CorrelationMatrix = CorrelationMatrix(1:obj.Nr.SamplesTotal^2,1:obj.Nr.SamplesTotal^2);
        end
        
        function [TimeCorrelation,Time] = PlotTimeCorrelation(obj,TimeSpacing)
            % Plot the time correlation. The input argument represents
            % the time-spacing and plots discrete points at multiples
            % of the time-spacing (only for presentational purpose)
            
            [TimeCorrelation,Time]=obj.GetTimeCorrelation;
            figure();
            plot(Time,abs(TimeCorrelation));
            ylabel('|Time Correlation|');
            xlabel('Time (s)');
            
            if not(exist('TimeSpacing','var'))
                TimePointSymbols = (-ceil(Time(end)/TimeSpacing):ceil(Time(end)/TimeSpacing))*TimeSpacing;
                hold on;
                Plot1 = stem(TimePointSymbols,abs(interp1(Time,TimeCorrelation,TimePointSymbols)),'black');
                legend(Plot1,{'TimeSpacing'});            
            end    
        end
        
        function [FrequencyCorrelation,Frequency] = PlotFrequencyCorrelation(obj,FrequencySpacing)
            % Plot the frequency correlation. The input argument represents
            % the time-spacing and plots discrete points at multiples
            % of the time-spacing (only for presentational purpose)
            
            [FrequencyCorrelation,Frequency]=obj.GetFrequencyCorrelation;
            figure();
            plot(Frequency,abs(FrequencyCorrelation));
            ylabel('|Frequency Correlation|');
            xlabel('Frequency (Hz)');
            
            if exist('FrequencySpacing','var')
                FrequencyPointSymbols = (-ceil(Frequency(end)/FrequencySpacing):ceil(Frequency(end)/FrequencySpacing))*FrequencySpacing;
                hold on;
                Plot1 = stem(FrequencyPointSymbols,abs(interp1(Frequency,FrequencyCorrelation,FrequencyPointSymbols)),'black');
                legend(Plot1,{'FrequencySpacing'});
            end
        end
    end
end

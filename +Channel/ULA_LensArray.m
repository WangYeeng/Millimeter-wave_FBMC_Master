function [effective_channel] = ULA_LensArray(N_BS,N_user,N_user_antenna,channel)
% This function is used to transform the spatial channel into beamspace channel
U_BS = zeros(N_BS);
for ii = 1:N_BS
    U_BS(:,ii) = (1/sqrt(N_BS))*Channel.steering_vec(N_BS,ii*1/N_BS);
end
effective_channel = [];
for kk = 1:N_user
    effective_channel_k = [];
    U_user = [];
    for ii = 1:N_user_antenna  % each user utilize a lens antenna
        U_user = [U_user 1/sqrt(N_user_antenna)*Channel.steering_vec(N_user_antenna,((ii-1)-(N_user_antenna-1)/2)*pi/N_user_antenna)]; % the resolution at the bs is N_BS
    end
    % generate effective channel
    effective_channel_k = U_user'*channel((kk-1)*N_user_antenna + 1:kk*N_user_antenna,:)*U_BS';
    effective_channel = [effective_channel; effective_channel_k];
end
return;
end


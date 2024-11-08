function [effective_channel] = UPA_LensArray(N_ver,N_hor,N_user,channel)
% this function is used to transform the spatial channel into beamspace channel
U_BS = zeros(N_ver*N_hor);
for ii = 1:N_hor
    for jj = 1:N_ver
        steering_aod_ver = 1/sqrt(N_ver)*Channel.steering_vec(N_ver,jj/N_ver);
        steering_aod_hor = 1/sqrt(N_hor)*Channel.steering_vec(N_hor,ii/N_hor);
        U_BS(:,((ii-1)*N_ver+jj)) = kron(steering_aod_hor,steering_aod_ver);
    end
end


effective_channel = [];
for kk = 1:N_user
    effective_channel_k = [];
    % generate effective channel
    effective_channel_k = channel(kk,:)*U_BS;
    effective_channel = [effective_channel; effective_channel_k];
end
return;
end


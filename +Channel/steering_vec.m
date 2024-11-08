function vector = steering_vec(N,theta)
%this function generates a steering vector according to the transmit
%antennas and angle of departure
vector = zeros(N,1);
for ii = 0:N-1
    vector(ii+1) = exp(1j*2*pi*theta*(ii - (N-1)/2));
end
end


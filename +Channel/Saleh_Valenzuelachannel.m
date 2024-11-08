function [h] = Saleh_Valenzuelachannel(N_ver,N_hor,NrPath,type)
         % Returns a Saleh_Valenzuelachannel matrix.
         % Supports only uniform linear and planar arrays (type = 'ULA' or 'UPA'). 

         M = N_ver*N_hor;
         d = 0.5;
         h = zeros(M,1);
         alpha = (normrnd(0, 1, NrPath, 1) + 1i*normrnd(0, 1, NrPath, 1)) / sqrt(2);

         switch type
               case 'ULA'
                   phi   = pi*rand(1,NrPath)-pi/2;
                   for l = 1:NrPath
                       a = 1/sqrt(M)*exp(-1i*2*pi*(-(M-1):2:M-1)'/2*d*sin(phi(l)));  
                       h = h + alpha(l)*a;
                   end
               case 'UPA'
                   phi1   = pi*rand(1,NrPath)-pi/2;
                   phi2   = pi*rand(1,NrPath)-pi/2;
                   for l  = 1:NrPath
                       a1 = 1/sqrt(N_ver)*exp(-1i*2*pi*(-(N_ver-1):2:N_ver-1)'/2*d*sin(phi1(l))*sin(phi2(l)));
                       a2 = 1/sqrt(N_hor)*exp(-1i*2*pi*(-(N_hor-1):2:N_hor-1)'/2*d*cos(phi2(l)));
                       a  = kron(a1,a2);
                       h  = h + alpha(l)*a;
                   end
         end
         h = sqrt(M/NrPath)*h;
end


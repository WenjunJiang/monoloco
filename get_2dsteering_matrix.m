function S = get_2dsteering_matrix(theta_vals,d_vals,n_sub,f,df,ant_sep,n_ant) 
%% SpotFi 2D Steering Matrix for the given theta_vals and d_vals
% Inputs:
% theta_vals: the search space for Angle of Arrival
% d_vals: the search space for Time of Flight
% n_sub: number of subcarriers
% f: Center Frequency of the transmission
% df: subcarrier bandwidth
% ant_sep: antenna sepration on the AP
% n_ant: number of receiving antennas
% Ouputs:
% S: 2D Steering vector for SpotFi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi = zeros(n_ant-1,length(theta_vals));
for i = 1:length(theta_vals)
    phi_theta = exp(-1j*2*pi*f/3e8*sin(theta_vals(i))*ant_sep);
    for j=1:n_ant-1
        Phi(j,i) = phi_theta.^(j-1);
    end
end
Omega = zeros(n_sub/2,length(d_vals));
for i = 1: length(d_vals)
    omega_t = exp(-1j*2*pi*df*d_vals(i)/3e8);
    Omega(:,i) = omega_t.^((1:n_sub/2)');
end

S = kron(Phi,Omega);
end
function [aoa,tof] = xy_to_aoa_tof(ap,xy)
%% This function converts the XY labels to AoA-ToF w.r.t AP
% Input:
% ap_loc : XY location of the AP antennas in meters N_antx2 location vector
% xy : XY locations of ground truth labels Nx2 locations for N location
% labels
% Ouputs:
% AoA : Angle of arrival of the N XY label points, NX1 vectors
% AoA : Time of flight of the N XY label points, NX1 vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aoa = zeros(n_points,1);
tof = zeros(n_points,1);
parfor i=1:n_points
    ap_pos = mean(ap);
    ap_vec=ap(1,:)-ap(end,:);
    X=xy(i,1)-ap_pos(1);
    Y=xy(i,2)-ap_pos(2);
    aoa(i)= rad2deg(sign(sum([X,Y].*ap_vec))*...
        (pi/2-acos(abs(sum([X,Y].*ap_vec))/norm([X,Y])/norm(ap_vec))));
    tof(i) = norm([X,Y]-ap_vec);
end
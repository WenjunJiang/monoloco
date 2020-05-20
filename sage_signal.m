function latent_signal = sage_signal(latent_parameter, latent_index, sage_const)
% SAGE_SIGNAL: Calculate latent multipath signal based on estimated parameters.
    T = sage_const.T;
    F = sage_const.F;
    A = sage_const.A;
    AoD = sage_const.D;
    
    tof_matrix = permute(repmat(sage_const.tof_candidates(:,latent_index(1)), 1, T, A, AoD), [2,1,3,4]);
    aoa_matrix = permute(repmat(sage_const.aoa_candidates(:,latent_index(2)), 1, T, F, AoD), [2,3,1,4]);
    aod_matrix = permute(repmat(sage_const.aod_candidates(:,latent_index(3)), 1, T, F, A), [2,3,4,1]);
    doppler_matrix = repmat(sage_const.doppler_candidates(:,latent_index(4)), 1, F, A, AoD);
    latent_signal = latent_parameter(5) .* tof_matrix .* aoa_matrix .* aod_matrix .* doppler_matrix;
end
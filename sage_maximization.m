function [latent_parameter, latent_index] = sage_maximization(latent_signal, latent_parameter, latent_index, sage_const)
% SAGE_MAXMIZATION: Maximization step of SAGE.
% LATENT_SIGNAL: T x F x A x AoD matrix
% CURRENT_PARAMETER: 1 x 5 vector: tof, aoa, doppler, amplitude
    
    T = size(latent_signal, 1);
    F = size(latent_signal, 2);
    A = size(latent_signal, 3);
    AoD = size(latent_signal, 4);
    doppler_range = sage_const.doppler_range;
    aoa_range = sage_const.aoa_range;
    aod_range = sage_const.aod_range;
    tof_range = sage_const.tof_range;

    % Estimation of tof
    aoa_matrix = permute(repmat(sage_const.aoa_candidates(:,latent_index(2)), 1, T, F, AoD), [2,3,1,4]);
    aod_matrix = permute(repmat(sage_const.aod_candidates(:,latent_index(3)),1,T,F,A),[2,3,4,1]);
    doppler_matrix = repmat(sage_const.doppler_candidates(:,latent_index(4)), 1, F, A, AoD);
    
    coeff_matrix = latent_signal .* conj(aoa_matrix) .* conj(aod_matrix) .* conj(doppler_matrix);
    coeff_vector = reshape(sum(sum(sum(coeff_matrix, 1), 4), 3), F, 1);
    coeff_vector = repmat(coeff_vector, 1, length(tof_range));
    object_vector = abs(sum(coeff_vector .* conj(sage_const.tof_candidates),1));
    [~, latent_index(1)] = max(object_vector);
    latent_parameter(1) = tof_range(latent_index(1));
    
    % Estimation of aoa
    tof_matrix = permute(repmat(sage_const.tof_candidates(:,latent_index(1)), 1, T, A, AoD), [2,1,3,4]);
    
    coeff_matrix = latent_signal .* conj(doppler_matrix) .* conj(tof_matrix) .* conj(aod_matrix);
    coeff_vector = reshape(sum(sum(sum(coeff_matrix, 1), 2), 4), A, 1);
    coeff_vector = repmat(coeff_vector, 1, length(aoa_range));
    object_vector = abs(sum(coeff_vector .* conj(sage_const.aoa_candidates), 1));
    [~, latent_index(2)]= max(object_vector);
    latent_parameter(2) = aoa_range(latent_index(2));
    
    
    % Estimation of aod
    aoa_matrix = permute(repmat(sage_const.aoa_candidates(:,latent_index(2)), 1, T, F, AoD), [2,3,1,4]);
    
    coeff_matrix = latent_signal .* conj(doppler_matrix) .* conj(tof_matrix) .* conj(aoa_matrix);
    coeff_vector = reshape(sum(sum(sum(coeff_matrix, 1), 2), 3), AoD, 1);
    coeff_vector = repmat(coeff_vector, 1, length(aod_range));
    object_vector = abs(sum(coeff_vector .* conj(sage_const.aod_candidates), 1));
    [~, latent_index(3)]= max(object_vector);
    latent_parameter(3) = aod_range(latent_index(3));
    
    % Estimation of doppler
    aod_matrix = permute(repmat(sage_const.aod_candidates(:,latent_index(3)),1,T,F,A),[2,3,4,1]);
    
    coeff_matrix = latent_signal .* conj(aoa_matrix) .* conj(tof_matrix) .* conj(aod_matrix);
    coeff_vector = reshape(sum(sum(sum(coeff_matrix, 2), 4), 3), T, 1);
    coeff_vector = repmat(coeff_vector, 1, length(doppler_range));
    object_vector = abs(sum(coeff_vector .* conj(sage_const.doppler_candidates), 1));
    [~, latent_index(4)] = max(object_vector);
    latent_parameter(4) = doppler_range(latent_index(4));
    
    % Estimation of amplitude
    doppler_matrix = repmat(sage_const.doppler_candidates(:,latent_index(4)), 1, F, A, AoD);
    coeff_matrix = latent_signal .* conj(aoa_matrix) .* conj(aod_matrix) .* conj(tof_matrix) .* conj(doppler_matrix);
    latent_parameter(5) = sum(coeff_matrix(:)) / (T * F * A * AoD);
end
function [final_parameter, residue_error] = sage_sfg(csi_signal, initial_parameter, initial_index, sage_const)
% SAGE_SFG: SAGE algorithm.
% CSI_SIGNAL: T x F x A matrix
    T = sage_const.T;
    F = sage_const.F;
    A = sage_const.A;
    AoD = sage_const.D;
    L = sage_const.L;
    N = sage_const.N;
    update_ratio = sage_const.update_ratio;
    
    % Initialize
    latent_signal = zeros(T, F, A, AoD, L);
    for ii = 1:L
        if initial_parameter(5,ii) ~= 0
            latent_signal(:,:,:,:,ii) = sage_signal(initial_parameter(:,ii), initial_index(:,ii), sage_const);
        end
    end
    
    % Iteration
    final_parameter = initial_parameter;
    temp_parameter = initial_parameter;
    temp_index = initial_index;
    final_index = initial_index;
    for ii = 1:N
        for jj = 1:L
            temp_signal = sage_expectation(csi_signal, latent_signal, jj, update_ratio);
            [temp_parameter(:,jj), temp_index(:,jj)] = sage_maximization(temp_signal, final_parameter(:,jj), final_index(:,jj), sage_const);
            latent_signal(:,:,:,:,jj) = sage_signal(temp_parameter(:,jj), temp_index(:,jj), sage_const);
        end
        parameter_diff = sqrt(sum(abs(temp_parameter - final_parameter).^2, 2));
        final_parameter = temp_parameter;
        final_index = temp_index;
        if parameter_diff(1) < 1e-9 && parameter_diff(2) < 1 / 180 * pi && parameter_diff(3) < 1 && parameter_diff(4) < 1e-9
            break;
        end
    end
    residue_error = csi_signal - sum(latent_signal, 4);
    residue_error = mean(abs(residue_error(:)))/mean(abs(csi_signal(:)));
end
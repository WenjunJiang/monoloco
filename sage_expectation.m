function expect_signal = sage_expectation(csi_signal, latent_signal, expect_index, update_ratio)
% SAGE_EXCEPTATION: Expectation step of SAGE.
% CSI_SIGNAL: T x F(=30) x A(=3) x AoD(=3) matrix
% LATENT_SIGNAL: T x F x A x AoD x L matrix
% EXPECT_INDEX: \in [1,L]
% UPDATE_RATIO: beta = 1
    noise_signal = csi_signal - sum(latent_signal, 5);
    expect_signal = latent_signal(:,:,:,:,expect_index) + update_ratio * noise_signal;
end
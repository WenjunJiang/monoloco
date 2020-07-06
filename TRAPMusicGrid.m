function [mus] = TRAPMusicGrid(aTot,csi_data,n_iter,delayLoS,aoaLoS,aodLoS,AS_vec,DS_vec,calib_const)
% % % Input
% aTot: Equivalent to L_scan in TRAP. We can treat nSubcarr*nRx as sensors
%       Since each position has a unique (tof,aoa), we can treat all
%       possible combination of tof and aoa as source topographies
% csi_data: (nSubcarr*nRx, nSample). We use it mainly to calculate covariance
%       matrix C=UDU' and signal space projection P=U_sU_s'
% n_iter: be 1 if we only want to remove LoS.
% % % Output
% mus: music spectrum, scanning-function value for all scanning space,
%       with dimension (n_iter,n_scan)
% % % RAP-music vs TRAP-music
% TRAP-music is different from RAP-music in Uso,Us,Uk, where the columns
% are truncated
L_scan = aTot;
[n_sens,n_scan] = size(L_scan);
C_meas = csi_data*csi_data';
mus = zeros(n_iter, n_scan);
B = zeros(n_sens, n_iter); % measurement of each source, equivalent to Ahat

remain_dim = 5;
%SVD & space of the measurement covariance matrix
[Utemp,~,~] = svd(C_meas);
Uso = Utemp(:,1:n_iter+remain_dim); %Here +1 because n_iter does not include LoS path
% Ahat: the measurement of LoS, equivalent to B(:,1) and l_found
Ahat = gridSampleBackscatter3D(delayLoS,aoaLoS,aodLoS,AS_vec,DS_vec,calib_const);
B(:,1) = Ahat;
% Calculated from LoS
for ITER = 1:n_iter
    l_found = B(:,1:ITER);
    %check for the conditioning
    [U,S] = svd(l_found'*l_found,'econ');
    s = diag(S);
    tol = 10*eps(s(1))*ITER; %a bit lower tolerance than in 'pinv'
    keep = sum(s>tol);
    if keep<ITER
        %truncated pseudoinverse and out-projection
        fprintf('trapscan_presetori: Ill conditioning in out-projection at iteration %d; truncating.\n',ITER);
        fprintf('                    You are likely trying to find more sources than the model supports.\n');
        Qk = eye(n_sens) - l_found*U(:,1:keep)*diag(1./s(1:keep))*U(:,1:keep)'*l_found';
    else
        %normal pseudoinverse and out-projection
        Qk = eye(n_sens) - l_found*U*diag(1./s)*U'*l_found';
    end
    
    %apply out-projection to forward model
    L_this = Qk*L_scan;
    [Us,~,~] = svd(Qk*Uso,'econ');
    %TRAP truncation
    Uk = Us(:,1:(n_iter-ITER+remain_dim));
    %scan over all test sources, matrix form.
    %in principle we would use a subspace projector Ps = Uk*Uk', but as
    %only norms are needed, it is faster to use Uk directly.
    L_this_normsq = sum(L_this.*L_this,1); % same to RAP-music
    PsL_normsq = sum((Uk'*L_this).^2,1); % same to RAP-music
    
    mus(ITER,:) = PsL_normsq./L_this_normsq;
    mus(ITER,:) = abs(mus(ITER,:)).^2;
end

end


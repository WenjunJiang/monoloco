function [steeringVec] = gridSampleBackscatter3D(delay,aoa,aod,AS_vec,DS_vec,calib_const)
% % % Input
% 
% % % Output
%
cc = calib_const;
aoaSteering = exp(-1j*2*pi*AS_vec'*cc.fc/cc.C*cos(aoa));
aodSteering = exp(-1j*2*pi*DS_vec'*cc.fc/cc.C*cos(aod));
delaySteering = exp(-1j*2*pi*cc.f'*delay);
steeringVecDelayAoATmp = delaySteering*aoaSteering.';
steeringVecDelayAoA = steeringVecDelayAoATmp(:);
steeringVecTmp = steeringVecDelayAoA*aodSteering.';
steeringVec = steeringVecTmp(:);

end


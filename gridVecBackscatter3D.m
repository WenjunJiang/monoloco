%% creates set of steering vectors. MUSIC spectrum is calculated at each of these steering vectors.
function [aTot] = gridVecBackscatter3D(delayRange,AoARange,AoDRange,AS_vec,DS_vec,calib_const,generateAtot)
% % % Input
% delayRange: possible values of delay
% AoARange: possible values of AoA
% AoDRange: possible values of AoD
% generateAtot: whether we need to generate a new Atot or load the previous one
% calib_const: parameters
% % % Output
% aTot: matrix with each column being a steering vector Tx_Rx_nSubCarr * DR_AR_AoDR
nTx = size(DS_vec,2); % number of antennas on Tx
nRx = size(AS_vec,2); % number of antennas on Rx

cc = calib_const;
nSubCarr = size(cc.f,2);

numGridPoints = length(delayRange)*length(AoARange)*length(AoDRange); % all combinations
if ~generateAtot
    load('aTotSave');
else
    aTot = zeros(nSubCarr*nRx*nTx, numGridPoints);
    i = 1;
    for aod = AoDRange
        for aoa = AoARange
            for delay = delayRange
                steeringVec = gridSampleBackscatter3D(delay,aoa,aod,AS_vec,DS_vec,calib_const);
                aTot(:,i) = steeringVec;
                i=i+1;
            end
        end
    end
    save('aTotSave','aTot')
end

end


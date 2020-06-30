function [music_spectrum] = RAPMusicGrid(aTot,Qn,Qs,delayLoS,aoaLoS,aodLoS,AS_vec,DS_vec,useNoise,calib_const)
if ~useNoise
    Ahat = gridSampleBackscatter3D(delayLoS,aoaLoS,aodLoS,AS_vec,DS_vec,calib_const);
    PerpAhat = eye(size(Ahat,1)) - Ahat*pinv(Ahat'*Ahat)*Ahat';
    Qs = PerpAhat*Qs;
    doHadamard = true;
    if doHadamard
        A = PerpAhat * aTot;
        QA = Qs'*A;
        music_spectrum = norms(QA,2,1)./norms(A,2,1);
        music_spectrum = music_spectrum.^2;
    end
end

end


%% obtain signal and noise subspace in MUSIC algorithm
function [Pn,Ps,Qn,Qs,EigenInfo] = GetQnBackscatter(X, nComps)
% % % Input
% X: 2 dimensions, nRx*nTx_nSubCarr*nSample
% nComps: number of paths
% % % Output
% Qn: matrix where each column is a basis vector of noise subspace (read MUSIC algorithm)
% Qs: matrix where each column is a basis vector of signal subspace (read MUSIC algorithm)                                               
% Pn = (Qn*Qn');
% Ps = (Qs*Qs');
% EigenInfo: structure containing some useful information like all the eigenvalues.

% % by singular value
% [U,D,~] = svd(X);
% % by eigenvalue
[Utmp,D] = eig(X*X');
D = abs(D);
[Dtmp,I] = sort(diag(D), 'descend');
D = diag(Dtmp);
U = Utmp(:,I);

minMP = 2;
useMDL = 1;
% % % MDL criterion based
MDL = [];
lambdaTot = diag(D);
subarraySize = size(X,1);
nSegments = size(X,2);
maxMultipath = length(lambdaTot); % previously used 6 for maximum number of multipath
for k = 1:maxMultipath
    MDL(k) = -nSegments*(subarraySize-(k-1))*log(geo_mean(lambdaTot(k:end))/mean(lambdaTot(k:end))) + 0.5*(k-1)*(2*subarraySize-k+1)*log(nSegments);
end
% % Another attempt to take the number of multipath as minimum of MDL
[~, SignalEndIdxTmp] = min(MDL);
if useMDL
    SignalEndIdx = max(SignalEndIdxTmp-1, 1);
end
SignalEndIdx = max(SignalEndIdx, minMP);

if ~isempty(nComps)
    SignalEndIdx = nComps; 
end

% % % % % Plotting the figures to demnstrate difference in eigenvalues.
% figure(1); plot(diff(db(diag(D))),'d-')
% % % SignalEndIdx = 10;
% figure(2); plot((diag(D)),'d-')
% figure(3); plot(zeros(length(diag(D)),1), diag(D),'d');
% hold on;
% % plot([-1 1], 0.01*[1 1], 'r--')
% plot([-1 1], 0.1*D(1)*[1 1], 'r--')
% hold off;
% figure(4); plot(MDL, 'd-');
% pause(0.1)
% % sprintf('min criterion is %d', find(Criterion1 & Criterion3,1,'last'))
% sprintf('Noise space dimesion is %d',size(U,2)-SignalEndIdx)
% sprintf('Signal space dimesion is %d',SignalEndIdx)
% sprintf('min MDL index is %d', SignalEndIdxTmp)
% 
Qn = U(:,SignalEndIdx+1:end);
Pn = (Qn*Qn');

Qs = U(:,1:SignalEndIdx);
Ps = (Qs*Qs');

EigenInfo = struct;
EigenInfo.Umatrix = U;
EigenInfo.nMPs = SignalEndIdx;
EigenInfo.singularValuesdB = db(diag(D));

end


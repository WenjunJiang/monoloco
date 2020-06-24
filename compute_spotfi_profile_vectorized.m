function P = compute_spotfi_profile_vectorized(h, theta_vals, d_vals, threshold,S)
%% Our Vectorized implementation SpotFi
% INPUT 
% h: 4 times 234 matrix of csi measurements (with slope across subcarriers
% removed)
% theta_vals: values of time where the profile has to be evaluate
% d_vals: values of distance where the profile has to be evaluated
% threshold: Threshold value to determine the noise space and signal space
% OUTPUT 
% p: is a length(theta_vals) times length(d_vals) array of complex
% values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Steps
% Create the Signal matrix
% Find eigen values
% Detect Noise Subspace
% Compute projects on the noise subspace
% disp('Entered vectorized spotfi');

h=h.';

[n_ant,n_sub] = size(h);

if(size(h,2)~=n_sub)
    fprintf('h does not have 30 subcarriers. Check the code\n');
end
A = zeros((n_ant-1)*n_sub/2, n_sub); %super resolution.
for i=1:n_sub/2
    for j=1:n_ant-1
        A((j-1)*n_sub/2 +(1:n_sub/2),i) = h(j,i-1+(1:n_sub/2)).';
    end
end
for i=n_sub/2+1:n_sub
    for j=2:n_ant
        A((j-2)*n_sub/2 +(1:n_sub/2),i) = h(j,i-n_sub/2-1+(1:n_sub/2)).';
    end
end
R = A*A';

[V, D] = eig(R);
eig_vals= diag(D);
idx = find(eig_vals<threshold*max(abs(eig_vals)));

As = 1./vecnorm((V(:,idx))'*S,2).^2; %This is the MUSIC spectrum

% Reshape the AoA-ToF Profile matrix
P = reshape(As,length(d_vals),length(theta_vals));
P = P.';

end
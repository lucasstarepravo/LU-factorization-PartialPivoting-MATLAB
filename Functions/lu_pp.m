function [A, P] = lu_pp(A)
n = length(A); P = 1:n;
for k = 1:n-1
    [~, max_index] = max(abs(A(k:n,k)));             % Gathering row of max index w.r.t. k
    max_index = max_index + k - 1;                   % Globalizing row index
    A([k max_index],:) = A([max_index k],:);         % Permuting A matrix
    P([k max_index])   = P([max_index k]);           % Permuting vector
    gauss_w = A(k+1:n,k)/A(k,k);                     % Computing Gaussian Weights
    A(k+1:n,k:n) = A(k+1:n,k:n) - gauss_w.*A(k,k:n); % Gaussian elimination
    A(k+1:n,k) = gauss_w';                           % Substituting GW in L part of A matrix
    end
end
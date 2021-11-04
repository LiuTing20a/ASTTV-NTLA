function [X] = prox_tensorLaplace(Z, rho, Xold, epsilon)
[n1,n2,n3] = size(Z);
n12 = min(n1,n2);

Z = fft(Z,[],3);
Xold = fft(Xold,[],3);

U = zeros(n1,n12,n3);
V = zeros(n2,n12,n3);
S = zeros(n12,n12,n3);

for i = 1 : n3
    [U(:,:,i),s,V(:,:,i)] = svd(Z(:,:,i),'econ');
    w = svd(Xold(:, :, i), 'econ');
    
    %laplace
    w = 1/epsilon*exp(-w/epsilon);
       
    s = diag(s);
    s = max(s-w*rho,0);
    %   s = max(s-rho*n3*w,0);
    S(:,:,i) = diag(s);
    
    
    X(:,:,i) = U(:,:,i)*S(:,:,i)*V(:,:,i)';
end
X = ifft(X,[],3);


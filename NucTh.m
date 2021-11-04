function X=NucTh(X,lambda,r)
if isnan(lambda)
    lambda=0;
end
[M,N,L] = size(X);
X = reshape(X,M*N,L);
[u,s,v]= svd(X,0);
s = diag(s);
if r<min(size(X))
    s(r+1:end) = 0;
end
s1 = max(s - lambda,0);
X=u*diag(s1)*v';
X = reshape(X,M,N,L);
end


function esti_rank = rank_estimation(X);
%
% The function returns an estimated rank of the input Matrix
%
% created by Fanhua Shang on 12/06/2016, fhshang@cse.cuhk.edu.hk
%
[m, n] = size(X);
num = 100;
if num>min(m,n)
    num = min(m,n);
end
S0 = svds(X, num);
S1 = S0(1:end-1)-S0(2:end);
S2 = S1./mean(S1(end-5:end));

% Spectral gap
r1 = 0; lam = 0.03; 
while(r1 <= 0)
    for i = 1:length(S2)
        ratio(i) = lam*max(S2(i:end)) + i;
    end
    [a1,idx1] = min(ratio);
    r1  = max(idx1-1);
    lam = lam + 0.03;
end
clear ratio a1 idx1;

if (r1<6)
    r11 = 0; lam = 0.03; 
    while (r11<=0)
        for i = 1:length(S2)-r1
            ratio(i) = lam*max(S2(i+r1:end))+i+r1;
        end
        [a1, idx1] = min(ratio);
        r11 = max(idx1-1)+r1;
        lam = lam + 0.03; 
    end
    r1 = r11;
end
clear ratio a1 idx1;

[a2,idx2] = max(S2);
if (idx2<5)
    [a2,idx2] = max(S2(5:end));  %%% rank >= 5
    r2 = idx2+4;
else
    r2 = idx2;
end
clear a2 idx2;

% Misssing data case 
if nnz(X)>m*n
    esti_rank = min(max(r1,r2),30);
else
    kappa = sqrt(m*n)/nnz(X);
    for i = 1:length(S0)-1
        ratio(i) = (S0(i+1) + S0(1)*sqrt(i*kappa))/S0(i);
    end
    [a3,idx3] = min(ratio);
    r3 = max(idx3);
    clear ratio a3 idx3;
    esti_rank = min(max([r1,r2,r3]),30);
end
end

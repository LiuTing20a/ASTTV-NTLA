function [X, his] = laplace_tsvd(B, Omega, opts)
tol = 1e-8; 
maxit = 500;
epsilon = 1e-5;
rho = 1.1;
beta = 1e-4;
max_beta = 1e10;
DEBUG = 0;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    maxit = opts.maxit;    end
if isfield(opts, 'epsilon');         epsilon = opts.epsilon;              end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'beta');          beta = opts.beta;                end
if isfield(opts, 'max_beta');      max_beta = opts.max_beta;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
if isfield(opts, 'Xtrue');      Xtrue = opts.Xtrue;         end

szB = size(B);
X = B; 
% X(logical(1-Omega)) = mean(B(Omega)); %%令初始估计中的丢失像素区域像素值平均化

% Y = zeros(szB); %% 辅助变量
Y = ones(szB);
% Y = rand(szB);
% Y = B;
% M = zeros(szB); %% 拉格朗日乘子
% M = ones(szB);
% M = rand(szB);
M = B;

his = [];

for k = 1:maxit
    Xold = X;
    %%% solve Y-subproblem
    [Y] = prox_tensorLaplace(Xold + M/beta, 1/beta, Y,epsilon);
    
    %%% solve X-subproblem
    X = Y - M/beta;
    
    X(Omega) = B(Omega);
    
    psnrlap(k) = myPSNR(Xtrue(:), X(:));
        
%     imshow(X), drawnow

    %%% check the convergence
    if isfield(opts, 'Xtrue')
        real(k) = norm(X(:) - Xtrue(:)) / norm(Xtrue(:));
    end
    res(k) = norm(X(:) - Xold(:)) / norm(Xold(:));
    if mod(k, 20) == 0 && DEBUG == 1
        fprintf('laplace: iter = %d   diff=%f\n', k, res(k));
    end
    if res(k) < tol && k > 40
        break;
    end
    
    %%% update Lagrange multiplier
    M = M + beta * (X-Y);  
    beta = min(rho * beta, max_beta);
end
his.res = res;
if isfield(opts, 'Xtrue')
    his.real = real;
end
his.iter = k;
his.psnrlap = psnrlap;
fprintf('laplace_tsvd ends: total iterations = %d   difference=%f\n', k, res(k));
end 





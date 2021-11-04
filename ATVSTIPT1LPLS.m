function [tenB, tenT,tenN, change] = ATVSTIPT1LPLS(tenD, lambda1,lambda2,lambda3,weight, mu,C,p,opts)

% Solve the Tensor Robust Principal Component Analysis based on weighted Tensor Nuclear Norm problem by ADMM
%
% min_{tenB,tenT} ||tenB||_Wb,*+lambda*||Wt.*tenT||_1, s.t. tenD=tenB+tenT
% initialize
epsilon         = 1e-6;
max_iter = 100;
normD = norm(tenD(:));
if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end


ext1 = 0;
extd = 0;
tenD = padarray(tenD,[ext1 ext1 extd],'symmetric');
[M1,N1,L1] = size(tenD);
lambda = 0.01;
pro =1.25;% tuning
lambda_s = 0.3;
AL_iters = 100;
mu=0.005;%0.01
r = 5;
ro = 1.2;  mu_max = 1e6;
TV_TYPE = 'iso';
im_size = [M1,N1,L1];
m1 = [1 0;0 0 ];
m2 = [-1 0 ;0 0 ];
template_time(:,:,1) = m1;
template_time(:,:,2) = m2;

FDx = psf2otf([1 -1],im_size); % fourier transform of Dx in 3D cube
FDy = psf2otf([1;-1],im_size); % fourier transform of Dx in 3D cube
FDz = psf2otf(template_time,im_size); % fourier transform of Dx in 3D cube

FDxH = conj(FDx);%用于计算复数的共轭值
FDyH = conj(FDy);
FDzH = conj(FDz);

% inverse
 IL= 1./((abs(FDx).^2 + abs(FDy).^2 + abs(FDz).^2) + 2);
%%  Initialization
sizeD= size(tenD);
B               = rand(sizeD);      % B : low-rank background image
Z               = B;                % Z : auxiliary variable for X
T               = zeros(sizeD);     % T : sparse target 
N               = zeros(sizeD);     % T : sparse target 
Y1              =zeros(sizeD);
%Lagrange Multipliers
d1 = zeros(M1,N1,L1); v1 = d1;
d2 = zeros(M1,N1,L1); v2 = d2;
d3 = zeros(M1,N1,L1); v3 = d3;
d4 = d3;
d5 = d3;
tol1 = 10e-8;
tol2  =10e-8;
%%
dim = size(tenD);
n=min(dim(1),dim(2));
change=zeros(1,max_iter);
for iter = 1 : max_iter
    Bk = B;
    Tk = T;
%% Update Z
    [Z] = prox_tensorLaplace(B-Y1/mu, 1/mu, Z,epsilon);
%% Update low-rank background tensor B and weightensor for B
     Lx =tenD-T-N+Z;
     Fx = IL.* (fftn(Lx+d1/mu) + FDxH.*fftn(v1 + d3/mu) + FDyH.*fftn(v2 + d4/mu) + FDzH.*fftn(v3 + d5/mu));
     B = ifftn(Fx);
    
%% Update sparse target tensor T and weightensor for T
    T = softthre(tenD-B+d1/mu-N,lambda2/mu);
%% updata V1,V2,V3
    v1 = MySoftTh( ifftn(FDx.*Fx) - d3./mu,lambda/mu );
    v2 = MySoftTh( ifftn(FDy.*Fx) - d4./mu,lambda/mu );
    v3 = MySoftTh( ifftn(FDz.*Fx) - d5./mu,lambda*pro/mu );
%% Update N
     N          = (mu*(tenD -B -T) + d1)/(mu+2*lambda3); 
%% Update multipliers Y1,Y2,Y3
    d1 = d1 + mu*(tenD -B - T-N);
    d2 = d2 + mu*(B - Z);
    d3 = d3 + mu*(v1- ifftn(FDx.*Fx));
    d4 = d4 + mu*(v2- ifftn(FDy.*Fx));
    d5 = d5 + mu*(v3- ifftn(FDz.*Fx));
    
%% Update mu
    mu = min(ro*mu,mu_max);
%% Stop Criterion 
    errList    = norm(tenD(:)-B(:)-T(:)-N(:)) / normD;
    fprintf('STTV: iterations = %d   difference=%f\n', iter, errList);
    change(iter)=(errList);
    if errList < epsilon
        break;  
    end 
end
%% Output
tenB=B;tenT=T;tenN=N;
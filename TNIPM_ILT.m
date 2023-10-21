function [x,status,history] = TNIPM_ILT(A,y,lambda)
%
% Multidimensional ILT l1-Regularized LS Problem Solver via TNIPM
% 
%
%   l1_ls solves problems of the following form:
%
%       minimize   ||A*x-y||^2 + lambda*sum(x_i),
%       subject to x_i >= 0, i=1,...,n
%
%   where A and y are problem data and x is variable (described below).
%
% INPUT
%   A       : the comprehensive Laplace Kernal Matrix, defined as the
%   Kroneck product of K1,K2,...,Kn
%   m       : number of  rows of A
%   n       : number of columns of A
%   y       : the vectorization of the acuqired data.
%   lambda  : positive scalar; regularization parameter
%
%
%   (advanced arguments)
%       eta     : scalar; parameter for PCG termination (default: 1e-3)
%       pcgmaxi : scalar; number of maximum PCG iterations (default: 5000)
%
% OUTPUT
%   x       : n vector; classifier
%   status  : string; 'Solved' or 'Failed'
%
%   history : matrix of history data. columns represent (truncated) Newton
%             iterations; rows represent the following:
%            - 1st row) gap
%            - 2nd row) primal objective
%            - 3rd row) dual objective
%            - 4th row) step size
%            - 5th row) pcg iterations
%            - 6th row) pcg status flag
%
% USAGE EXAMPLES
%   [x,status] = l1_ls_nonneg(A,y,lambda);
%   [x,status] = l1_ls_nonneg(A,At,m,n,y,lambda,0.001);
%
 
% Orignal AUTHOR:    Kwangmoo Koh <deneb1@stanford.edu>
% Reference: 
% Kim, S. J.; Koh, K.; Lustig, M.; Gorinevsky, S. B. D., An Interior-Point Method for Large-Scale L1-Regularized Least Squares. IEEE Journal of Selected Topics in Signal Processing 2008, 1 (4), 606-617.
% COPYRIGHT 2008 Kwangmoo Koh, Seung-Jean Kim, and Stephen Boyd


% Some Simplification was made by Enping Lin for the specific ILT inverse problem.
% UPDATE    January 10 2021
% If you use this code, please cite the related papers:
% Enping Lin, Ville-Veikko Telkki, Xiaoqing Lin, Chengda Huang, Haolin Zhan, Yu Yang*, Yuqing Huang*, and Zhong Chen*High-Resolution Reconstruction for Multidimensional Laplace NMR, J. Phys. Chem. Lett. 2021, 12, 21, 5085¨C5090.  https://doi.org/10.1021/acs.jpclett.1c01022


%------------------------------------------------------------
%       INITIALIZE
%------------------------------------------------------------

% IPM PARAMETERS
MU              = 2;        % updating parameter of t
MAX_NT_ITER     = 400;      % maximum IPM (Newton) iteration

% LINE SEARCH PARAMETERS
ALPHA           = 0.01;     % minimum fraction of decrease in the objective
BETA            = 0.5;      % stepsize decrease factor
MAX_LS_ITER     = 100;      % maximum backtracking line search iteration

% 
   At = A';
   [m,n] = size(A);

% VARIABLE ARGUMENT HANDLING
t0         = min(max(1,1/lambda),n/1e-3);
defaults   = {1e-8,true,1e-3,5000,ones(n,1),t0};
[reltol,quiet,eta,pcgmaxi,x,t] = deal(defaults{:});

f = -x;

% RESULT/HISTORY VARIABLES
pobjs = [] ; dobjs = [] ; sts = [] ; pitrs = []; pflgs = [];
pobj  = Inf; dobj  =-Inf; s   = Inf; pitr  = 0 ; pflg  = 0 ;

ntiter  = 0; lsiter  = 0; zntiter = 0; zlsiter = 0;
normg   = 0; prelres = 0; dx =  zeros(n,1);

% diagxtx = diag(At*A);
diagxtx = 2*ones(n,1);


%------------------------------------------------------------
%               MAIN LOOP
%------------------------------------------------------------

for ntiter = 0:MAX_NT_ITER
  
    z = A*x-y;
    
    %------------------------------------------------------------
    %       CALCULATE DUALITY GAP
    %------------------------------------------------------------

    nu = 2*z;

    minAnu = min(At*nu);
    if (minAnu < -lambda)
        nu = nu*lambda/(-minAnu);
    end
    pobj  =  z'*z+lambda*sum(x,1);
    dobj  =  max(-0.25*nu'*nu-nu'*y,dobj);%£¿
    gap   =  pobj - dobj;

    pobjs = [pobjs pobj]; dobjs = [dobjs dobj]; sts = [sts s];
    pflgs = [pflgs pflg]; pitrs = [pitrs pitr];

    %------------------------------------------------------------
    %   STOPPING CRITERION
    %------------------------------------------------------------

    if (gap/abs(dobj) < reltol) 
        status  = 'Solved';
        history = [pobjs-dobjs; pobjs; dobjs; sts; pitrs; pflgs];
        if (~quiet) disp('Absolute tolerance reached.'); end
          disp(sprintf('total pcg iters = %d\n',sum(pitrs)));
          return;
    end
    %------------------------------------------------------------
    %       UPDATE t
    %------------------------------------------------------------
    if (s >= 0.5)
        t = max(min(n*MU/gap, MU*t), t);
    end

    %------------------------------------------------------------
    %       CALCULATE NEWTON STEP
    %------------------------------------------------------------
    
    d1 = (1/t)./(x.^2);

    % calculate gradient
    gradphi = [At*(z*2)+lambda-(1/t)./x];
    
    % calculate vectors to be used in the preconditioner
    prb     = diagxtx+d1;

    % set pcg tolerance (relative)
    normg   = norm(gradphi);
    pcgtol  = min(1e-1,eta*gap/min(1,normg));%?
    
    if (ntiter ~= 0 && pitr == 0) pcgtol = pcgtol*0.1; end

if 1
    [dx,pflg,prelres,pitr,presvec] = ...
        pcg(@AXfunc_l1_ls,-gradphi,pcgtol,pcgmaxi,@Mfunc_l1_ls,...
            [],dx,A,At,d1,1./prb);%?
end
    %dx = (2*A'*A+diag(d1))\(-gradphi);

    if (pflg == 1) pitr = pcgmaxi; end
    
    %------------------------------------------------------------
    %   BACKTRACKING LINE SEARCH
    %------------------------------------------------------------
    phi = z'*z+lambda*sum(x)-sum(log(-f))/t;
    s = 1.0;
    gdx = gradphi'*dx;
    for lsiter = 1:MAX_LS_ITER
        newx = x+s*dx;
        newf = -newx;
        if (max(newf) < 0)
            newz   =  A*newx-y;
            newphi =  newz'*newz+lambda*sum(newx)-sum(log(-newf))/t;
            if (newphi-phi <= ALPHA*s*gdx)
                break;
            end
        end
        s = BETA*s;
    end
    if (lsiter == MAX_LS_ITER)
         break; 
    end % exit by BLS
        
    x = newx; f = newf;
end


%------------------------------------------------------------
%       ABNORMAL TERMINATION (FALL THROUGH)
%------------------------------------------------------------
if (lsiter == MAX_LS_ITER)
    % failed in backtracking linesearch.
    if (~quiet) disp('MAX_LS_ITER exceeded in BLS'); end
    status = 'Failed';
elseif (ntiter == MAX_NT_ITER)
    % fail to find the solution within MAX_NT_ITER
    if (~quiet) disp('MAX_NT_ITER exceeded.'); end
    status = 'Failed';
end
history = [pobjs-dobjs; pobjs; dobjs; sts; pitrs; pflgs];

 return;

%------------------------------------------------------------
%       COMPUTE AX (PCG)
%------------------------------------------------------------
function [y] = AXfunc_l1_ls(x,A,At,d1,p1)

y = (At*((A*x)*2))+d1.*x;

%------------------------------------------------------------
%       COMPUTE P^{-1}X (PCG)
%------------------------------------------------------------
function [y] = Mfunc_l1_ls(x,A,At,d1,p1)

y = p1.*x;

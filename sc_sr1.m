function [sigmaStar,sigma_perp,sigma_par,pStar,opt1,opt2,spd_check,N_iter] ...
    = sc_sr1( g, delta, gamma, Psi,invM);

%
% Copyright (2017): Johannes Brust, Jennifer Erway, Roummel Marcia
%
%
% The technical report and software are available at 
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
%
% This code is based on the manuscript entitled, "Algorithm xxx: SC-SR1:
%  MATLAB Software for Solving Shape-Changing L-SR1 Trust-Region
%  Subproblems"
%  by
% Johannes Brust, Oleg Burdakov, Jennifer Erway, Roummel Marcia, and
% Ya-xiang Yuan
%
% This code is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%---------------------------------------------------------------------------    
% 
% Inputs: 
%         g: the gradient vector and right-hand side of the first
%                optimality condition
%         Psi, invM - used to define the quasi-Newton matrix
%         delta: the trust-region radius
%         gamma: scaling of initial Hessian matrix approximation B0, 
%                i.e., B0 = gamma*I
%
% Output: sigmaStar: the optimal Lagrange multiplier for the 2-norm
%                       subproblem.
%         sigma_perp: The optimal Lagrange multiplier with the shape-changing norm (for the problem in v_perp).
%         sigma_par: The optimal Lagrange multiplier with the shape-changing norm (for the problem in v_par).
%         pStar: the optimal solution to the trust-region subproblem
%         opt1: the residual of the 1st optimality condition: norm(BpStar+sigmaStar*pStar + g) 
%         opt2: the residual of the 2nd optimality condition: sigmaStar*abs(delta-norm(pStar))
%         spd_check: gives the smallest eig of (B+sigmaStar I): Lambda(1) + sigmaStar 
%         N_iter: number of Newton iterations
%
% Note: the variables "show" and "altOpt" can be toggled as follows:
%       show    = 0: Runs silently and assigns no value to opt1, opt2, spd_check
%               = 1: Runs silently but assigns value to opt1, opt2, spd_check
%               = 2: Verbose setting; also assigns value to opt1, opt2, spd_check
%       altOpt  = 1: Computes the 2nd optimality condition as the sum of
%                       squares for the shape-changing norms.
%               = 2: Computes the 2nd optimality condition as the sum of
%                       absolute values for the shape-changing norms.

show        = 1;     % verbosity flag
altOpt      = 2;     % optimality condition 2 flag. 


%% Initializations
maxIter     = 100;   % maximum number of iterations for Newton's method
tol         = 1e-10; % tolerance for Newton's method
N_iter      = 0;     % number of Newton iterations

% Initializing optimal Lagrange multipliers.
sigmaStar   = 0; 
sigma_perp  = 0; 
sigma_par   = 0; 


%% Ensure Psi is full rank
if (rank(Psi)~=size(Psi,2))
  fprintf('\n\n');
  ME = MException('Catastropic_error:Psi',' Psi is not full rank... exiting obs.m');
  throw(ME);
end

%% Compute eigenvalues Hat{Lambda} using Choleksy
PsiPsi      = Psi'*Psi;
R           = chol(PsiPsi);
RMR         = R* (invM\R'); 
RMR         = (RMR+(RMR)')./2;  %symmetrize RMR', if needed
[U D ]      = eig(RMR);

%% Eliminate complex roundoff error then sort eigenvalues and eigenvectors
[D,ind]     = sort(real(diag(D)));  %D is now a vector
U           = U(:,ind);            

% Compute Lambda as in Equation (7) and Lambda(1)
n           = size(g,1);
sizeD       = size(D,1);
Lambda       = D + gamma.*ones(sizeD,1);  %formerly, Lambda_one
Lambda      = Lambda.*(abs(Lambda)>tol); %thresholds

% Define P_parallel and g_parallel 
RU          = R\U;     
P_parallel  = Psi*RU;  
Psig        = Psi'*g;  
g_parallel  = RU'*Psig; 

% Compute a_j = (g_parallel)_j for j=1...k+1; a_{k+2}=||g_perp||
a_j         = g_parallel;
a_kp2      = sqrt(abs(g'*g - g_parallel'*g_parallel));

if a_kp2 < tol  %fix small values of a_kp2 
  a_kp2 = 0;
end

if Lambda(1)<0   
  rp0 = find((Lambda-Lambda(1))>tol,1);
  vmin_parallel = -a_j(rp0:end)./(Lambda(rp0:end)-Lambda(1));
else
  rp1 = find(Lambda>0,1);
  v0_parallel = a_j(rp1:end)./Lambda(rp1:end);
end

% Compute sigma_perp
if (abs(a_kp2/gamma) <= delta) && (gamma > 0)
  sigma_perp  = 0;
else
  sigma_perp  = a_kp2/delta - gamma;
  if abs(a_kp2) < sqrt(tol) % ||g_perp|| = 0
    sigma_perp = -gamma;
  end
end

% Compute beta, where v_perp* = beta*g_perp.
% Compute pStarPerp, where pStar = pStarPar + pStarPerp
beta        = 0;
pStarPerp   = zeros(n,1);

if sigma_perp + gamma ~= 0

  if abs(a_kp2) > sqrt(tol)
    beta        = -1/(sigma_perp+gamma);
    pStarPerp   = beta*(g + P_parallel*(-g_parallel));
  end
else
  beta    = delta;
  e       = zeros(n,1);
  found   = 0; 
  for i=1:sizeD
    e(i)=1;
    u_min = e + P_parallel*(-P_parallel(i,1:end))';
    if norm(u_min)>tol
      found = 1;
      break;
    end
    e(i)=0;
  end  
  if found == 0
    e(m+1) = 1;
    u_min = e + P_parallel*(-P_parallel(i,1:end))';
  end
  u_min       = u_min/(norm(u_min));
  pStarPerp   = beta*u_min;
end

if ((Lambda(1)>0) & (norm(v0_parallel)<=delta))
  
  sigmaStar = 0;

elseif (Lambda(1)<=0) & (phiBar_f(-Lambda(1),delta,Lambda,a_j)>0)

  sigmaStar = -Lambda(1);
  sigma_par = -Lambda(1); % Extension with shape-changing norm
  
  % forms v = (Lambda_one + sigmaStar I)^\dagger P_\parallel^Tg
  index_pseudo      = find(abs(Lambda+sigmaStar)>tol);
  v                 = zeros(sizeD+1,1);
  v(index_pseudo)   = a_j(index_pseudo)./(Lambda(index_pseudo)+sigmaStar); 
  
  % forms pStar using Equation (16)
  % Shape-changing extension
  pStarPar    = P_parallel*(-v(1:sizeD));
  
  if (abs(gamma+sigmaStar)<tol) 
    pStar = pStarPar;
  else
    pStar = pStarPar + pStarPerp;
  end
  
  if Lambda(1) < 0
    
    alpha = sqrt(delta^2-pStarPar'*pStarPar);
    pHatStar = pStar;
    
    % compute z* using Equation (17)
    if abs(Lambda(1)-Lambda(1))<tol 
      zstar = (1/norm(P_parallel(:,1))).*alpha*P_parallel(:,1); 
    else %gamma=Lambda(1)
      e     = zeros(size(g,1),1);
      found = 0; 
      for i=1:sizeD
	e(i)    =1;
	u_min   = e+ P_parallel*(-P_parallel(i,1:end))';
	if norm(u_min)>tol
	  found = 1;
	  break;
	end
	e(i)=0;
      end  
      if found == 0
	e(m+1)  = 1;
	u_min   = e+ P_parallel*(-P_parallel(i,1:end))';
      end
      u_min = u_min/(norm(u_min));
      zstar = alpha*u_min;  
    end
    
    pStar = pHatStar+zstar;
    
  end  
else
  if Lambda(1)>0
    [sigmaStar, N_iter] = Newton(0,maxIter,tol,delta,Lambda,a_j);
  else
    sigmaHat = max(abs(a_j)/delta - Lambda);
    if sigmaHat>-Lambda(1) 
      [sigmaStar, N_iter] = Newton(sigmaHat,maxIter,tol,delta,Lambda,a_j);
    else
      [sigmaStar, N_iter] = Newton(-Lambda(1),maxIter,tol,delta,Lambda,a_j);
    end
  end
  
  %Shape-changing extension
  sigma_par = sigmaStar; 
  pStar     = P_parallel*(-a_j./(Lambda+sigmaStar)) + pStarPerp;

end


%%% optimality check
%% Checks extended for the Shape-Changing norms
if show>=1
  
  PpStar      = P_parallel'*pStar;
  BpStar      = (gamma+sigma_perp).*pStar + Psi*(invM\(Psi'*pStar)) + (sigma_par-sigma_perp).*P_parallel*(PpStar);
  
  opt1        = norm(BpStar + g);
  
  nPstarsq    = pStar'*pStar;
  nPpStarsq   = PpStar'*PpStar;
  
  opt21       = sigma_par*(sqrt(nPpStarsq)-delta);
  opt22       = sigma_perp*(sqrt(nPstarsq-nPpStarsq)-delta);
  opt2        = opt21 + opt22; 
  
  if altOpt == 1
    opt2 = opt21^2 + opt22^2;
  elseif altOpt == 2
    opt2 = abs(opt21)+ abs(opt22);
  end
  
  spd_check       = zeros(2,1);
  spd_check(1)    = Lambda(1) + sigma_par;
  spd_check(2)    = gamma + sigma_perp;
  
  
  if show==2
    fprintf('\nOptimality condition #1: %8.3e', opt1);
    fprintf('\nOptimality condition #2: %8.3e', opt2);

      fprintf('\nLambda(1)+sigma_par: %8.2e', spd_check(1));
      fprintf('\ngamma+sigma_perp: %8.2e', spd_check(2));
    
    fprintf('\n\n');
  end
else
  opt1 = [ ];
  opt2 = [ ];
  spd_check = [ ];
      spd_check = zeros(2,1);

  phiBar_check = [ ];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ phiBar ] = phiBar_f( sigma,delta,D,a_j )
%phiBar_f evaluates the continous extension of 
%phi(sigma) = 1/ ||v|| - 1/delta. 
%
%For further details, please see the following technical report:
%
% "OBS: MATLAB Solver for L-SR1 Trust-Region Subproblems"
% by Johannes Burst, Jennifer Erway, and Roummel Marcia
%
% Copyright (2015): Johannes Brust, Jennifer Erway, and Roummel Marcia
%
% The technical report and software are available at 
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
%
% This code is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%---------------------------------------------------------------------------    


m = size(a_j,1); 
D = D + sigma.*ones(m,1);   %vector
eps_tol = 1e-10;

% test if numerator or denominator has a zero
if ( sum( abs(a_j) < eps_tol ) > 0 ) | (sum( abs(D) < eps_tol ) > 0 )    
  pnorm2 = 0;
  for i = 1:m        
    if (abs(a_j(i)) > eps_tol) & (abs(D(i)) < eps_tol)
      phiBar = -1/delta;
      return; 
    elseif (abs(a_j(i)) > eps_tol) & (abs(D(i)) > eps_tol)
      pnorm2    = pnorm2 + (a_j(i)/D(i))^2;
    end
  end
  phiBar     = sqrt(1/pnorm2) - 1/delta;
  return;
end

%% numerators and denominators are nonzero
p = a_j./D;
phiBar = 1/norm(p) - 1/delta;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ phiBar, phiBar_g ] = phiBar_fg( sigma,delta,D,a_j )
%phiBar_f evaluates the continous extension of 
%phi(sigma) = 1/ ||v|| - 1/delta and its derivative.
%
%
% Copyright (2015): Johannes Brust, Jennifer Erway, and Roummel Marcia
%
% The technical report and software are available at 
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
%
% This code is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%---------------------------------------------------------------------------    

m        = size(a_j,1); 
D        = D + sigma.*ones(m,1);   %vector
eps_tol  = 1e-10; 
phiBar_g = 0;

% test if numerator or denominator has a zero
if ( sum( abs(a_j) < eps_tol ) > 0 ) || (sum( abs((D)) < eps_tol ) > 0 )    
  pnorm2 = 0;
  for i = 1:m        
    if (abs(a_j(i)) > eps_tol) && (abs(D(i)) < eps_tol)
      phiBar   = -1/delta; 
      phiBar_g = 1/eps_tol;
      return; 
    elseif abs(a_j(i)) > eps_tol && abs(D(i)) > eps_tol
      pnorm2   = pnorm2   +  (a_j(i)/D(i))^2;
      phiBar_g = phiBar_g + ((a_j(i))^2)/((D(i))^3);
    end
  end
  normP    = sqrt(pnorm2);
  phiBar   = 1/normP - 1/delta;
  phiBar_g = phiBar_g/(normP^3);
  return;
end

%%% Numerators and denominators are all nonzero
% Compute phiBar(sigma)
p      = a_j./D;
normP  = norm(p);
phiBar = 1/normP - 1/delta;

phiBar_g = sum((a_j.^2)./(D.^3));
phiBar_g = phiBar_g/(normP^3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ x,k ] = Newton(x0,maxIter,tol,delta,Lambda,a_j)
%Newton finds a zero of phiBar.
%
%For further details, please see the following technical report:
%
% Copyright (2015): Johannes Brust, Jennifer Erway, and Roummel Marcia
%
% The technical report and software are available at 
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
%
% This code is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%---------------------------------------------------------------------------    


x = x0;  %initialization
k = 0;   %counter

[f g]  = phiBar_fg(x,delta,Lambda,a_j);
f0     = eps*norm(f);  %relative error tolerance
tau    = 1e-8;  %absolute error tolerance
N_flag = 0;


while (N_flag==0) & (k < maxIter)
    x       = x - f / g;
    [f g]  = phiBar_fg(x,delta,Lambda,a_j);
    k = k + 1;
    if (norm(f) <= (norm(f0)+tau))
      N_flag = 1;
    end
end
if k>50
 fprintf('\ntoo many iterations.. pausing');
 keyboard
end



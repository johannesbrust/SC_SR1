function [ ] = driver_sc_sr1(S,Y,gamma,g,delta,varargin)
% driver for sc_sr1.m
% to run with a file: driver_sc_sr1([ ],[ ],[ ],[ ],[ ],'filename.mat');
% to run without a file: driver_sc_sr1(S,Y,gamma,g,delta);
%
% Input: S - change in the iterates
%        Y - change in the gradients of the iterates
%        gamma - initialization B0=gamma I
%        g  - the righthand side
%        delta - the trust-region radius
%
%        optional: a file in single quotes that includes Psi, invM, gamma, g, delta.
%
%        One must either supply S, Y, gamma, g, delta or the optional
%        file name
%
% Copyright (2017): Johannes Brust, Jennifer Erway, Roummel Marcia
%
% This code is based on the manuscript entitled, "Algorithm xxx: SC-SR1:
%  MATLAB Software for Solving Shape-Changing L-SR1 Trust-Region
%  Subproblems"
%  by
% Johannes Brust, Oleg Burdakov, Jennifer Erway, Roummel Marcia, and
% Ya-xiang Yuan
%
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


if nargin==6
  experiment_data = varargin{1};
  load(experiment_data);   %load data
  n = size(Psi,1);
else
  n = size(S,1);
  SY          = S'*Y;
  SS          = S'*S;
  invM        = tril(SY)+ tril(SY,-1)'-gamma.*SS;
  invM        = (invM+(invM)')./2;  %symmetrize it
  Psi         = Y-gamma.*S;
end
  
% initialize variables
kp1 = 5;  %number of limited memory updates
n_runs = 5;  % number of test runs for testing time  %5 
i_runs = 3;  % the run whose time is stored (i_runs<=n_runs) %3
tol  = 1e-10;

% do n_runs runs, save time associated with i_runs

for testruns = 1:n_runs
  tic
  [sigmaStar,sigma_perp,sigma_par,pStar,opt1,opt2,spd_check,N_iter]=sc_sr1(g,delta,gamma,Psi,invM);
  timer = toc;
  if testruns==i_runs
    save_time = timer;
  end
end

%%%output
fprintf('\nRepeating each problem %d times saving the time of run #%d...\n\n',n_runs,i_runs);
fprintf('\n   n       opt1      opt2    lambda1+sig_par gamma+sig_perp');
fprintf('   sig_par   sig_perp    N_iter ');
fprintf('  time');
fprintf('\n-------  --------  --------  --------------- --------------- --------- ');
fprintf(' ---------  -------  --------  ');

fprintf('\n%1.1e  %8.2e  %8.2e     ',n,opt1,opt2);
fprintf(' %8.2e    ', spd_check(1));
fprintf('    %8.2e ', spd_check(2));
fprintf('    %8.2e  ', sigma_par);
fprintf(' %8.2e ', sigma_perp);

fprintf(' %5d    ', N_iter);
fprintf('%8.2e  \n ', save_time);   
fprintf('\n\n');



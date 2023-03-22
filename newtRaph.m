%% Simple Implementation of the Newton-Raphson-Method with constraints
%
% Mapping of x in R^n --> fun in R^m
%
% (c) Michael Mauersberger 2021 (v0.1), 2023 (v1.0), LGPL v2.1
%
% Newton-Raphson method with constraints. Damping coefficient helps to find
% a feasible solution. Number of sweep points help to set a new initial
% searching point on the axis between the minimum and maximum of the
% domain. Search ends if function tolerance or argument tolerance has been
% reached. Linear constraints are applied via writing back the arguments to
% the constraint (done by means of QR decomposition and a pseudo-inverse).
%
% References:
% Wikipedia - "Newtonverfahren" (German), "Newton's method" (English)
%

function [xVal,fVal,iter,exit] = newtRaph(fun,x0_in,xLow_in,xUpp_in,  ...
    ALim_in,bLim_in,fTol_in,xTol_in,maxIter_in,dampCoeff_in,methDiff_in,hDiff_in,nSweep_in)

% lower bound
xLow = [];
% upper bound
xUpp = [];
% linearly defined constraints
ALim = [];
bLim = [];
% function tolerance value
fTol = 1e-8;
% argument tolerance value
xTol = eps;
% maximum iteration number
maxIter = 100;
% damping coefficient
% --> enhances solution behaviour, cf. example problem:
% [xVal,fVal,it,ex] = newtRaph(@(x)(x.^3-2*x+2),0,[],[],[],[],1e-10,1e-10,100,.0,'cent',1e-10,2)
% [xVal,fVal,it,ex] = newtRaph(@(x)(x.^3-2*x+2),0,[],[],[],[],1e-10,1e-10,100,.2,'cent',1e-10,2)
dampCoeff = .2;
% method of calculating the difference quotient (forw, back, cent)
methDiff = 'forw';
% step size of difference quotient
hDiff = 1e-8;
% number of sweep points on out-of-limit restore
nSweep = 2;

% number of correction iterations
nCorr = 10;

if nargin > 1
  if ~isempty(x0_in)
    x0 = x0_in(:);
  end
end
if nargin > 2
  if ~isempty(xLow_in)
    xLow = xLow_in(:);
  end
end
if nargin > 3
  if ~isempty(xUpp_in)
    xUpp = xUpp_in(:);
  end
end
if nargin > 4
  if ~isempty(ALim_in)
    ALim = ALim_in;
  end
end
if nargin > 5
  if ~isempty(bLim_in)
    bLim = bLim_in(:);
  end
end
if nargin > 6
  if ~isempty(fTol_in)
    fTol = fTol_in;
  end
end
if nargin > 7
  if ~isempty(xTol_in)
    xTol = xTol_in;
  end
end
if nargin > 8
  if ~isempty(maxIter_in)
    maxIter = maxIter_in;
  end
end
if nargin > 9
  if ~isempty(dampCoeff_in)
    dampCoeff = dampCoeff_in;
  end
end
if nargin > 10
  if ~isempty(methDiff_in)
    methDiff = methDiff_in;
  end
end
if nargin > 11
  if ~isempty(hDiff_in)
    hDiff = hDiff_in;
  end
end
if nargin > 12
  if ~isempty(nSweep_in)
    nSweep = nSweep_in;
  end
end

% number of arguments
nx = length(x0);

% check input parameters
try
  fun(x0);
catch
  error('Function form and argument mismatch!')
end
if ~isempty(xLow) && nx ~= length(xLow) || ~isempty(xUpp) && nx ~= length(xUpp)
  error('Limits form mismatch!')
end
if xor(isempty(ALim),isempty(bLim))
  error('Linear constraints incomplete!')
elseif ~isempty(ALim)
  if size(ALim,2) ~= nx || size(ALim,1) ~= length(bLim)
    error('Linear constraints form mismatch!')
  end
end
if isnumeric(fTol)
  fTol = fTol(1);
else
  error('Wrong type for fTol!')
end
if isnumeric(xTol)
  xTol = xTol(1);
else
  error('Wrong type for xTol!')
end
if isnumeric(maxIter)
  maxIter = fix(maxIter(1));
else
  error('Wrong type for maxIter!')
end
if isnumeric(dampCoeff)
  dampCoeff = dampCoeff(1);
else
  error('Wrong type for dampCoeff!')
end
if ~ischar(methDiff)
  error('Wrong type for methDiff!')
end
if isnumeric(hDiff)
  hDiff = hDiff(1);
else
  error('Wrong type for hDiff!')
end
if isnumeric(nSweep)
  nSweep = fix(nSweep(1));
else
  error('Wrong type for nSweep!')
end
% minimum number of sweep points is 2 (lower and upper bound)
nSweep = max(2,nSweep);

% initial function evaluation
xVal = x0;
xVal_old = [];
fVal = fun(xVal);

% start periodically at the opposite border if no solution has been found
bordSweep = [0 0];
iter = 0;
exit = 0;

while true
  
  % check domain (lower and upper limits)
  lBndTest = true(nx,1);
  uBndTest = true(nx,1);
  if ~isempty(xLow) && ~isempty(xUpp)
    if ~isempty(xLow)
      lBndTest = (xVal >= xLow);
    end
    if ~isempty(xUpp)
      uBndTest = (xVal <= xUpp);
    end
    if all(~lBndTest) || all(~uBndTest)
      if sum(bordSweep) == nSweep
        warning('No roots subject to boundaries has been found!')
        exit = -2;
        break
      end
      if all(~lBndTest)
        bordSweep(1) = bordSweep(1) + 1;
        % ratio of sweep between minimum and maximum point in domain
        wSweep = (bordSweep(1) - 1) / (nSweep - 1);
        xVal = (1 - wSweep) * xLow + wSweep * xUpp;
      end
      if all(~uBndTest)
        bordSweep(2) = bordSweep(2) + 1;
        % ratio of sweep between minimum and maximum point in domain
        wSweep = (bordSweep(2) - 1) / (nSweep - 1);
        xVal = (1 - wSweep) * xUpp + wSweep * xLow;
      end
    else
      if any(~lBndTest)
        xVal(~lBndTest) = xLow(~lBndTest);
      end
      if any(~uBndTest)
        xVal(~uBndTest) = xUpp(~uBndTest);
      end
    end
  end

  % check linear constraints
  if ~isempty(ALim) && ~isempty(bLim)
    linLimTest = all(ALim * xVal <= bLim);
    if ~linLimTest
      % QR decomposition
      [~,RLim,iCols] = qr(ALim);
      iiCols = arrayfun(@(x)(find(iCols(x,:))),(1:size(iCols,1)));
      iZero = all(abs(RLim) < eps,1);
      indepR = find(~iZero);
      % find write-back value via pseudo-inverse
      xValinv = xVal + pinv(ALim) * (bLim - ALim * xVal);
      xVal(iiCols(indepR)) = xValinv(iiCols(indepR));
    end
  end
  
  % function evaluation with updated x
  fVal = fun(xVal);
  fVal = fVal(:);
  
  % if function/argument tolerance or maximum number of iterations has been reached
  % argument tolerance only if values are inside domain
  bxTol = all(lBndTest) && all(uBndTest) && (~isempty(xVal_old) && norm(xVal - xVal_old) <= xTol);
  if norm(fVal) <= fTol || bxTol || iter >= maxIter
    if iter >= maxIter
      warning('Maximum iteration number of %d reached!',maxIter)
      exit = 1;
    end
    if norm(fVal) > fTol
      if norm(xVal - xVal_old) <= xTol
        warning('Argument change lower than tolerance!')
      end
      warning('No roots within tolerance has been found!')
      exit = -1;
    end
    break
  end
  
  % correction with hDiff if Jacobian is zero
  dx_corr = hDiff;
  while true
    J = jacobian(fun,xVal,hDiff,methDiff);
    if any(J(:) ~= 0)
      break
    else
      warning('Jacobian is zero! Adding %g',dx_corr)
      % small variation by means of random correction value
      xVal = xVal + dx_corr;
      dx_corr = dx_corr * 2;
      if dx_corr > hDiff * 2^nCorr
        warning('No roots found due to zero Jacobian!')
        exit = -2;
        break
      end
    end
  end
  if exit < 0
    break
  end
  % solve new step
  if rank(J) < max(size(J))
    dx = pinv(J) * -fVal;
  else
    dx = J \ -fVal;
  end
  
  xVal_old = xVal;
  xVal = xVal_old + (1 - dampCoeff) * dx;
  iter = iter + 1;
  
end

% column vector
xVal = xVal(:);



% calculating the jacobian
function J = jacobian(funLoc,x,h,methDiff)

  f = funLoc(x);
  hlim = abs(f) * eps;
  if h < hlim && hlim < 1 / eps
    warning('Difference step corrected: %g --> %g',h,hlim)
    hcorr = hlim;
  else
    hcorr = h;
  end
  dimX = length(x);
  dimF = length(f);
  J = zeros(dimF,dimX);
  
  for i = 1:dimF
    for j = 1:dimX
      hdim = zeros(size(x));
      hdim(j) = 1;
      switch methDiff
        case 'cent'
          f1 = funLoc(x + hcorr * hdim);
          f2 = funLoc(x - hcorr * hdim);
          dh = 2 * hcorr;
        case 'forw'
          f1 = funLoc(x + hcorr * hdim);
          f2 = f;
          dh = hcorr;
        case 'back'
          f1 = f;
          f2 = funLoc(x - hcorr * hdim);
          dh = hcorr;
        otherwise
          error('Calculation method "%s" of difference quotient unknown!',methDiff)
      end
      if any(isnan(f1(:))) || any(isnan(f2(:))) || any(isinf(f1(:))) || any(isinf(f2(:)))
        % one as standard gradient
        J(i,j) = 1;
      else
        J(i,j) = (f1 - f2) / dh;
      end
    end
  end
  
end


end
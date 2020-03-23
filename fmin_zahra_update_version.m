function [x,fval,exitflag,output] = fmin_zahra_update_version(funfcn,x,lb,ub,options,varargin)
%FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%   X = FMINSEARCH(FUN,X0) starts at X0 and attempts to find a local minimizer 
%   X of the function FUN.  FUN is a function handle.  FUN accepts input X and 
%   returns a scalar function value F evaluated at X. X0 can be a scalar, vector 
%   or matrix.
%
%   X = FMINSEARCH(FUN,X0,OPTIONS)  minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created
%   with the OPTIMSET function.  See OPTIMSET for details.  FMINSEARCH uses
%   these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, FunValCheck,
%   PlotFcns, and OutputFcn.
%
%   X = FMINSEARCH(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminsearch' in PROBLEM.solver. 
%
%   [X,FVAL]= FMINSEARCH(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINSEARCH(...) returns an EXITFLAG that describes
%   the exit condition. Possible values of EXITFLAG and the corresponding
%   exit conditions are
%
%    1  Maximum coordinate difference between current best point and other
%       points in simplex is less than or equal to TolX, and corresponding 
%       difference in function values is less than or equal to TolFun.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSEARCH(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fminsearch(@sin,3)
%     finds a minimum of the SIN function near 3.
%     In this case, SIN is a function that returns a scalar function value
%     SIN evaluated at X.
%
%     FUN can be an anonymous function:
%        X = fminsearch(@(x) norm(x),[1;2;3])
%     returns a point near the minimizer [0;0;0].
%
%     FUN can be a parameterized function. Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
%        c = 1.5;                        % The parameter.
%        X = fminsearch(@(x) f(x,c),[0.3;1])
%        
%   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
%
%   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.

%   Copyright 1984-2017 The MathWorks, Inc.


defaultopt = struct('Display','notify','MaxIter','200*numberOfVariables',...
    'MaxFunEvals','200*numberOfVariables','TolX',1e-4,'TolFun',1e-4, ...
    'FunValCheck','off','OutputFcn',[],'PlotFcns',[]);

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && strcmpi(funfcn,'defaults')
    x = defaultopt;
    return
end

if nargin<3, options = []; end

% Detect problem structure input
if nargin == 1
    if isa(funfcn,'struct') 
        [funfcn,x,options] = separateOptimStruct(funfcn);
    else % Single input and non-structure
        error('MATLAB:fminsearch:InputArg',...
            getString(message('MATLAB:optimfun:fminsearch:InputArg')));
    end
end

if nargin == 0
    error('MATLAB:fminsearch:NotEnoughInputs',...
        getString(message('MATLAB:optimfun:fminsearch:NotEnoughInputs')));
end


% Check for non-double inputs
if ~isa(x,'double')
  error('MATLAB:fminsearch:NonDoubleInput',...
    getString(message('MATLAB:optimfun:fminsearch:NonDoubleInput')));
end

n = numel(x);
numberOfVariables = n;
%x = (lb + ub)/2;

% Check that options is a struct
if ~isempty(options) && ~isa(options,'struct')
    error('MATLAB:fminsearch:ArgNotStruct',...
        getString(message('MATLAB:optimfun:commonMessages:ArgNotStruct', 3)));
end

printtype = optimget(options,'Display',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
maxfun = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxiter = optimget(options,'MaxIter',defaultopt,'fast');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');

% In case the defaults were gathered from calling: optimset('fminsearch'):
if ischar(maxfun) || isstring(maxfun)
    if strcmpi(maxfun,'200*numberofvariables')
        maxfun = 200*numberOfVariables;
    else
        error('MATLAB:fminsearch:OptMaxFunEvalsNotInteger',...
            getString(message('MATLAB:optimfun:fminsearch:OptMaxFunEvalsNotInteger')));
    end
end
if ischar(maxiter) || isstring(maxiter)
    if strcmpi(maxiter,'200*numberofvariables')
        maxiter = 200*numberOfVariables;
    else
        error('MATLAB:fminsearch:OptMaxIterNotInteger',...
            getString(message('MATLAB:optimfun:fminsearch:OptMaxIterNotInteger')));
    end
end

switch printtype
    case {'notify','notify-detailed'}
        prnt = 1;
    case {'none','off'}
        prnt = 0;
    case {'iter','iter-detailed'}
        prnt = 3;
    case {'final','final-detailed'}
        prnt = 2;
    case 'simplex'
        prnt = 4;
    otherwise
        prnt = 1;
end
% Handle the output
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to plotfcns; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

header = ' Iteration   Func-count     min f(x)         Procedure';

% Convert to function handle as needed.
funfcn = fcnchk(funfcn,length(varargin));
% Add a wrapper function to check for Inf/NaN/complex values
if funValCheck
    % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
    % having to change the calls that look like this:
    % f = funfcn(x,varargin{:});
    % x is the first argument to CHECKFUN, then the user's function,
    % then the elements of varargin. To accomplish this we need to add the 
    % user's function to the beginning of varargin, and change funfcn to be
    % CHECKFUN.
    varargin = [{funfcn}, varargin];
    funfcn = @checkfun;
end

n = numel(x);

% Initialize parameters
% rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
% onesn = ones(1,n);
% two2np1 = 2:n+1;
% one2n = 1:n;

xM =x; % Middle point in the mentioned range.
fM = funfcn(xM,varargin{:});
changing = 2;
[xM,fM, xL, fL, xR, fR] = random_points(funfcn, xM, fM, lb, ub);


% Set up a simplex near the initial guess.
% xin = x(:); % Force xin to be a column vector
% v = zeros(n,n+1); fv = zeros(1,n+1);
v(:,1) = x;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
% x(:) = xin;    % Change x to the form expected by funfcn
fv(:,1) = fM;
func_evals = 1;
itercount = 0;
how = '';
% Initial simplex setup continues later

% Initialize the output and plot functions.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'init',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end

% Print out initial f(x) as 0th iteration
if prnt == 3
    disp(' ')
    disp(header)
    fprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how);
elseif prnt == 4
    formatsave.format = get(0,'format');
    formatsave.formatspacing = get(0,'formatspacing');
    % reset format when done
    oc1 = onCleanup(@()set(0,'format',formatsave.format));
    oc2 = onCleanup(@()set(0,'formatspacing',formatsave.formatspacing));
    format compact
    format short e
    disp(' ')
    disp(how)
    disp('v = ')
    disp(v)
    disp('fv = ')
    disp(fv)
    disp('func_evals = ')
    disp(func_evals)
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end

% Continue setting up the initial simplex.
% Following improvement suggested by L.Pfeffer at Stanford
% usual_delta = 0.05;             % 5 percent deltas for non-zero terms
% zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
% for j = 1:n
%     y = xin;
%     if y(j) ~= 0
%         y(j) = (1 + usual_delta)*y(j);
%     else
%         y(j) = zero_term_delta;
%     end
%     v(:,j+1) = y;
%     x(:) = y; f = funfcn(x,varargin{:});
%     fv(1,j+1) = f;
% end

% sort so v(1,:) has the lowest function value
% fv(:,1) = fM;
% fv(:,2) = fR;
% fv(:,3) = fL;
% v(:,1) = xM;
% v(:,2) = xR;
% v(:,3) = xL;
% 
% [fv,j] = sort(fv(1,:));
% v = v(:,j);

how = 'initial simplex';
itercount = itercount + 1;
func_evals = 1;
if prnt == 3
    fprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how)
elseif prnt == 4
    disp(' ')
    disp(how)
    disp('v = ')
    disp(v)
    disp('fv = ')
    disp(fv)
    disp('func_evals = ')
    disp(func_evals)
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end
exitflag = 1;

% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v3-v1||,...,||v(n+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the 
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations 
% are exceeded

fsaving_min = fv(:,1);
xsaving_min = v(:,1);
func_evals = func_evals+2;
    
while func_evals < maxfun && itercount < maxiter && fsaving_min > 0

    
%    fsaving_min = fM;
%    xsaving_min = xM;
    
    if (func_evals > maxfun)
        changing = 0;
        break
    end

    if (xL == xM & xR == xM)
        changing = 0;
        break
    end
    
%     if (changing == 0)
%         [xM,fM, xL, fL, xR, fR] = random_points(funfcn, xM, fM, lb, ub);
%     end
    
%     xNew1 = (xL+xM)/2;
%     xNew2 = (xM+xR)/2;
    
    %% Calculating Obj Fun for two new points:
%     fxNew1 = funfcn(xNew1);
%     fxNew2 = funfcn(xNew2);
    changing = 0;
    
    
    %% points have a V-shape; (xM,fM) is best (line 21)
    if (fM <= fL && fM <= fR)
        
        if (xM-xL) <= (xR-xM)
            
            xNew2 = (xR+xM)/2;
            fxNew2 = funfcn(xNew2);

            %% replacing xR with xNew2 if fM <= fxNew2
            if fM <= fxNew2
       
%                 fv(:,end) = fxNew2;
%             v(:,end) = xNew2;
%             [fv,j] = sort(fv(1,:));
%             v = v(:,j);
            
            if fM < fsaving_min
                
                fR = fxNew2;
                xR = xNew2;
                func_evals = func_evals+1;
                changing = 0;
            
                fsaving_min = fM;
                xsaving_min = xM;
                break;
            
                
                else
                    
                fR = fxNew2;
                xR = xNew2;
                func_evals = func_evals+1;
                changing = 1;
                
                end
                %% replacing xL with fM and fM with fxNew2
            else
                
%                 fv(:,end) = fxNew2;
%             v(:,end) = xNew2;
%             [fv,j] = sort(fv(1,:));
%             v = v(:,j);
            
            if fxNew2 < fsaving_min
                    changing = 0;
                    fsaving_min = fxNew2;
                    xsaving_min = xNew2; 
%                     fL = fM;
%                     xL = xM;
%                     fM = fxNew2;
%                     xM = xNew2;
                    func_evals = func_evals+1;

                    break
                    
                else
                    
                fL = fM;
                xL = xM;
                fM = fxNew2;
                xM = xNew2;
                func_evals = func_evals+1;
                changing = 1;
                
                end
            end
            
            %% if (xM-xL) <= (xR-xM) does not hold, then
        else
            xNew1 = (xL+xM)/2;
            fxNew1 = funfcn(xNew1);
            %% replacing xR with xNew1 if fM <= fxNew1
            if fxNew1 <= fM
                
%             fv(:,end) = fxNew1;
%             v(:,end) = xNew1;
%             [fv,j] = sort(fv(1,:));
%             v = v(:,j);
            
            if fxNew1 < fsaving_min
                    changing = 0;
                    fsaving_min = fxNew1;
                    xsaving_min = xNew1;
           
%                     fR = fM;
%                     xR = xM;
%                     fM = fxNew1;
%                     xM = xNew1;
                    func_evals = func_evals+1;
%                    changing = 1;
                    break
                   
                else
                    
                fR = fM;
                xR = xM;
                fM = fxNew1;
                xM = xNew1;
                func_evals = func_evals+1;
                changing = 1;
                
                end
                
            else
                
%                 fv(:,end) = fxNew1;
%             v(:,end) = xNew1;
%             [fv,j] = sort(fv(1,:));
%             v = v(:,j);
            
            if fM < fsaving_min
                    changing = 0;
                    fsaving_min = fM;
                    xsaving_min = xM;
%                 fL = fxNew1;
%                 xL = xNew1;
                func_evals = func_evals+1;
                %changing = 0;
               break
               
                else
                    
                fL = fxNew1;
                xL = xNew1;
                func_evals = func_evals+1;
                changing = 1;
                   
                end
            end
        end
        
        %% points have a /-, \-, or /\-shape; (xM,fM) is not best (line 27).
    elseif fL < fM
        
        xNew1 = (xL+xM)/2;
        fxNew1 = funfcn(xNew1);
        
        if fxNew1 < fL
            %% replacing fR with fM and fM with fxNew1
            fR = fM;
            xR = xM;
            fM = fxNew1;
            xM = xNew1;
            func_evals = func_evals+1;
            changing = 0;
            
            
            fv(:,1) = fM;
fv(:,2) = fR;
fv(:,3) = fL;
v(:,1) = xM;
v(:,2) = xR;
v(:,3) = xL;

[fv,j] = sort(fv(1,:));
v = v(:,j);
            
            if fv(:,1) < fsaving_min
                    
                    fsaving_min = fv(:,1);
                    xsaving_min = v(:,1);
            end
            changing = 0;
                    break;
            %% Else do nothing;
            
        else
            changing = 1;
            func_evals = func_evals+1;
            fR = fM;
            xR = xM;
            fM = fxNew1;
            xM = xNew1;
            %changing = 1;
%             if fM < fsaving_min
%             fsaving_min = fM;
%             xsaving_min = xM;
%             end
            % break;
            
            end
        
    else 
        xNew2 = (xR+xM)/2;
        fxNew2 = funfcn(xNew2);
        %% replacing fL with fM and fM with fxNew2
        if fxNew2 < fR
            
            fL = fM;
            xL = xM;
            fM = fxNew2;
            xM = xNew2;
            func_evals = func_evals+1;
            changing = 1;
            fv(:,1) = fM;
fv(:,2) = fR;
fv(:,3) = fL;
v(:,1) = xM;
v(:,2) = xR;
v(:,3) = xL;

[fv,j] = sort(fv(1,:));
v = v(:,j);
            changing = 0;
            if fv(:,1) < fsaving_min
                    
                    fsaving_min = fv(:,1);
                    xsaving_min = v(:,1);
            end
            break;
            %% Else do nothing;
        else
            fL = fM;
            xL = xM;
            fM = fxNew2;
            xM = xNew2;
            func_evals = func_evals+1;
            changing = 1;
%             if fM < fsaving_min
%             fsaving_min = fM;
%             xsaving_min = xM;
%             end
            %func_evals = func_evals+1;
            %changing = 0;
            %break;
            
            end            
            
        
        
        
    end
    fv = fsaving_min;
    v = xsaving_min;
%    [fv,j] = sort(fv);
%    v = v(:,j);
    itercount = itercount + 1;
    if prnt == 3
        fprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how)
    elseif prnt == 4
        disp(' ')
        disp(how)
        disp('v = ')
        disp(v)
        disp('fv = ')
        disp(fv)
        disp('func_evals = ')
        disp(func_evals)
    end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
            func_evals, how, fv(:,1),varargin{:});
        if stop  % Stop per user request.
            [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  prnt > 0
                disp(output.message)
            end
            return;
        end
    end
end   % while


%fsaving_min = fM;
%xsaving_min = xM;
x(:) = xsaving_min;
fval = fsaving_min;

output.iterations = itercount;
output.funcCount = func_evals;
output.algorithm = 'Zahra_Koen';

% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'done',itercount, func_evals, how, fval, varargin{:});
end

if func_evals >= maxfun
    msg = getString(message('MATLAB:optimfun:fminsearch:ExitingMaxFunctionEvals', sprintf('%f',fval)));
    if prnt > 0
        disp(' ')
        disp(msg)
    end
    exitflag = 0;
elseif itercount >= maxiter
    msg = getString(message('MATLAB:optimfun:fminsearch:ExitingMaxIterations', sprintf('%f',fval)));
    if prnt > 0
        disp(' ')
        disp(msg)
    end
    exitflag = 0;
else
    msg = ...
      getString(message('MATLAB:optimfun:fminsearch:OptimizationTerminatedXSatisfiesCriteria', ...
               sprintf('%e',tolx), sprintf('%e',tolf)));
    if prnt > 1
        disp(' ')
        disp(msg)
    end
    exitflag = 1;
end

output.message = msg;

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,state,iter,...
    numf,how,f,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.
%
% state - can have the values 'init','iter', or 'done'.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.iteration = iter;
optimValues.funccount = numf;
optimValues.fval = f;
optimValues.procedure = how;

xOutputfcn(:) = x;  % Set x to have user expected size
stop = false;
state = char(state);
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fminsearch:InvalidState',...
                getString(message('MATLAB:optimfun:fminsearch:InvalidState')));
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fminsearch:InvalidState',...
                getString(message('MATLAB:optimfun:fminsearch:InvalidState')));
    end
end

%--------------------------------------------------------------------------
function [x,FVAL,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT updates or sets all the output arguments of FMINBND when the optimization
% is interrupted.

% Call plot function driver to finalize the plot function figure window. If
% no plot functions have been specified or the plot function figure no
% longer exists, this call just returns.
callAllOptimPlotFcns('cleanuponstopsignal');

x = xOutputfcn;
FVAL = optimValues.fval;
EXITFLAG = -1;
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.algorithm = 'Nelder-Mead simplex direct search';
OUTPUT.message = getString(message('MATLAB:optimfun:fminsearch:OptimizationTerminatedPrematurelyByUser'));

%--------------------------------------------------------------------------
function f = checkfun(x,userfcn,varargin)
% CHECKFUN checks for complex or NaN results from userfcn.

f = userfcn(x,varargin{:});
% Note: we do not check for Inf as FMINSEARCH handles it naturally.
if isnan(f)
    error('MATLAB:fminsearch:checkfun:NaNFval',...
        getString(message('MATLAB:optimfun:fminsearch:checkfun:NaNFval', localChar( userfcn ))));  
elseif ~isreal(f)
    error('MATLAB:fminsearch:checkfun:ComplexFval',...
        getString(message('MATLAB:optimfun:fminsearch:checkfun:ComplexFval', localChar( userfcn ))));  
end

%--------------------------------------------------------------------------
function strfcn = localChar(fcn)
% Convert the fcn to a character array for printing

if ischar(fcn)
    strfcn = fcn;
elseif isstring(fcn) || isa(fcn,'inline')
    strfcn = char(fcn);
elseif isa(fcn,'function_handle')
    strfcn = func2str(fcn);
else
    try
        strfcn = char(fcn);
    catch
        strfcn = getString(message('MATLAB:optimfun:fminsearch:NameNotPrintable'));
    end
end



%% Generating two other points:
function [xM,fM, xL, fL, xR, fR] = random_points(funfcn, xM, fM, lb, ub)
%% Slope:
slope_up = 1; % lower range for slope. slope should be in Range [-1,1].
slope_down = -1; % upper range for slope. slope be in Range [-1,1].
slope_random = (slope_up - slope_down).*rand(1,10)+slope_down; % 10 random samples in Range [-1,1].
pos_slope_rand = randi(length(slope_random)); % Index (or position) of 10 random samples.
slope = slope_random(pos_slope_rand); % Pick one slope (according to the index is chosen )

%% Note: It is not necessary the random point starts passing M at the begining:
y = slope.* xM + 5 * rand(1, 50); % The linear line formula passing from middle point: y = slope * x + c. Generating 10 sample points.
%y(y < lb) = lb(y < lb);  y(y > ub) = ub(y > ub); % Force the sample points to be in the range.
pos_y = randi(length(y(1,:))); % Index (or position) of 10 random samples.
xR = y(:,pos_y); % Pick one point that pass M and xR.
%xR(xR < lb) = lb(xR < lb);  xR(xR > ub) = ub(xR > ub); % Force the sample points to be in the range.

%xL = xM - xR;
% y1 = slope.* xM + 5 * rand(1, 50); % The linear line formula passing from middle point: y = slope * x + c. Generating 10 sample points.
% %y(y < lb) = lb(y < lb);  y(y > ub) = ub(y > ub); % Force the sample points to be in the range.
% pos_y1 = randi(length(y(1,:))); % Index (or position) of 10 random samples.
% xL = y1(:,pos_y1); % Pick one point that pass M and xR.
% %xL(xL < lb) = lb(xL < lb);  xL(xL > ub) = ub(xL > ub); % Force the sample points to be in the range.


%% This part is written to map xR to one edge.
% uplow = {lb, ub}; % Take lower and upper line
% rand_uplow = randi([1, 2], 1); % Get a 1 or 2 randomly (place of them).
% thispoint = uplow(rand_uplow); % Get lower or upper.
% randomEdge = cell2mat(thispoint(:,1)); % Because the thispoint is saved in cell, transform it to mat.
% pos_low_up = randi(length(randomEdge(:,1))); % Index (or position) of one row of lower or upper (It depends which one is choden.)
% 
% if randomEdge == lb
%     xR(pos_low_up,1)= lb(pos_low_up,1); % If the lower is chosen, map xR to that.
% else
%     xR(pos_low_up,1)= ub(pos_low_up,1); % else the upper is chosen, map xR to that.
% end
% 
% % Finding the third point in the line that pass from M and xR:
dis = xM - xR; % the distance between M to xR.
xL = xM + dis; % New point is generated.

%% Calculating Obj Func for all three points:
fR = funfcn(xR);
fL = funfcn(xL);

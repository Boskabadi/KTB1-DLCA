function [sys,x0,str,ts,simStateCompliance] = RP_chrom(t,x,u,flag,gridsize)

switch flag,
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes(gridsize);
    %%%%%%%%%%%%%%%
    % Derivatives %
    %%%%%%%%%%%%%%%
    case 1,
        sys = mdlDerivatives(t,x,u,gridsize);
    %%%%%%%%%%
    % Update %
    %%%%%%%%%%
    case 2,
        sys = mdlUpdate(t,x,u);
    %%%%%%%%%%%%
    % Outputs % 
    %%%%%%%%%%%%
    case 3,
        sys = mdlOutputs(t,x,u,gridsize);
    %%%%%%%%%%%%%%%%%%%%%%%
    % GetTimeOfNextVarHit %
    %%%%%%%%%%%%%%%%%%%%%%%
    case 4,
        sys = mdlGetTimeOfNextVarHit(t,x,u);
    %%%%%%%%%%%%%
    % Terminate %
    %%%%%%%%%%%%%
    case 9,
        sys = mdlTerminate(t,x,u);
    %%%%%%%%%%%%%%%%%%%%
    % Unexpected flags %
    %%%%%%%%%%%%%%%%%%%%
    otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance] = mdlInitializeSizes(gridsize)

par = struct(...
    'components',3,...          % Number of components
    'nconc',3);                 % Number of concentrations

sizes = simsizes;
sizes.NumContStates  = par.nconc * par.components * gridsize;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = par.nconc * par.components * gridsize;
sizes.NumInputs      = par.components;
sizes.DirFeedthrough = 1;   % to understand whether inputs affect outputs (1=yes, 0=no)
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

% initialize the initial conditions
x0  = zeros(par.nconc * par.components * gridsize ,1); 
str = []; % str is always an empty matrix
ts  = [0 0]; % Initialize sample time

% Specify the block simStateCompliance
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys = mdlDerivatives(t,x,u,gridsize)

%% Parameters
par = struct(...
    'components',3,...          % Number of components
    'length', 500E-2,...        % Length of column [dm]
    'velocity', 20/1000/0.03*10,... % Linear velocity [dm/h] 
    'dax', 0.0113,...           % Axial dispersion [dm^2/h] 
    'epsbed',  0.35,...         % Bed porosity
    'epspar',  0.8,...          % Particle porosity
    'keff',  [0.4172 0.5606 0.3443],...  % Effective film transfer [dm/h] 
    'rad', 10E-6,...            % Particle radius [dm]
    'ads', [1E-07*0.9,1E-07*0.9,6*1.1E2]*3600E-6,... % Langmuir adsorption coefficient [dm^3/g/h]
    'des', [4E09*1.1,4E09*1.1,1*0.9E-3]*3600,... % Langmuir desorption coefficient [1/h]
    'qmax', [1E-6*0.9,1E-6*0.9,70*1.1]);      % Langmuir maximum concentration [g/dm^3]

nconc = 3; % Number of concentrations
h = par.length / (gridsize - 1); % Length of each cell

% Prepare sparse matrices only once
M = h/6 * diag(ones(gridsize-1,1),1) + ...
    2*h/3 * diag(ones(gridsize,1)) + ...
    h/6 * diag(ones(gridsize-1,1),-1);    
M(gridsize,gridsize) = h/3; 
M(1,1) = 1; 
M(1,2) = 0;
M = sparse(M);

C = (1/2) * diag(ones(gridsize-1,1), 1) + ...
    - (1/2) * diag(ones(gridsize-1,1),-1); 
C(1,1) = -1/2; 
C(gridsize,gridsize) =  1/2; 
C = sparse(C);

A = -1/h * diag(ones(gridsize-1,1), 1) + ...
    2/h * diag(ones(gridsize,1)) + ...
    -1/h * diag(ones(gridsize-1,1),-1);     
A(1,1) = 1/h; 
A(gridsize,gridsize) = 1/h; 
A = sparse(A);

% Reshape states and inputs for vectorized calculations
c = reshape(x, gridsize, par.components, nconc);
conc = u; % Directly use inputs from the tank

% Calculate derivatives using vectorized operations
dcdt = zeros(size(c)); % Preallocate derivative array
for j = 1:par.components
    % Mass balance for the first concentration
    dcdt(:,j,1) = -(par.velocity/par.epsbed * C + par.dax * A) * c(:,j,1) - ...
                  (1-par.epsbed)/par.epsbed * par.keff(j) * 3.0 / par.rad * M * (c(:,j,1) - c(:,j,2));

    % Boundary condition
    dcdt(1,j,1) = dcdt(1,j,1) - par.velocity/par.epsbed * (c(1,j,1) - conc(j));
    
    % Langmuir model for the third concentration
    sum = 1 - (c(:,3,1) / par.qmax(j)); % Simplified
    dcdt(:,j,3) = par.ads(j) * par.qmax(j) * sum .* c(:,j,2) - ...
                  par.des(j) * c(:,j,3);
    
    % Pore concentration
    dcdt(:,j,2) = -(1-par.epspar)/par.epspar * dcdt(:,j,3) + ...
                   par.keff(j) * 3 / par.rad / par.epspar * (c(:,j,1) - c(:,j,2));     
end

sys = reshape(dcdt, [], 1, 1); % Output as a column vector

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys = mdlUpdate(t,x,u)

sys = []; % No updates needed

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys = mdlOutputs(t,x,u,gridsize)

ncomp = 3; % Number of components
nconc = 3; % Number of concentrations

% Reshape the state vector to ensure it matches expected output size
sys = reshape(x, gridsize, ncomp, nconc); % Shape it based on gridsize

% Ensure the output is a real vector of length 450
if numel(sys) ~= 450
    error('Output must be a real vector of length 450. Current length: %d', numel(sys));
end

sys = sys(:); % Ensure output is a column vector of length 450

% end mdlOutputs


%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.
%=============================================================================
%
function sys = mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1; % Set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys = mdlTerminate(t,x,u)

sys = []; % No termination tasks needed

% end mdlTerminate

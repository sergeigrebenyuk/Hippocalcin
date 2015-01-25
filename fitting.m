clear all;
clear functions;

global fig_handle;

global simD1 fittedD1   % Diffusion coefficient (simulated and fitted) for cytosol
global simD2 fittedD2   % Diffusion coefficient (simulated and fitted) for plasma membrane
global simK1 fittedK1   % Membrane binding rate constant
global simK2 fittedK2   % Membrane unbinding rate constant
global Surf1  Surf2               % Simulated 2D-solutions 
global fittedSurf1 fittedSurf2    % Fitted 2D-solutions 
global simU1  fittedU1 rawU1 fakedU1 % 1D-solutions for cytosol at x=0: simulated fitted original and noised respetively
global simU2  fittedU2 rawU2 fakedU2
global noise1  noise2   % Noise added to simulated diffusion to produce 'real-looking' data
global roiSize;         % Size of FRAP'ed ROI (microns)
global modelLength      % Length of 1D dendrite (microns)
global deltaX           % Size of data grid cell (microns)
global simulationTime   % Simulation time (sec)
global deltaT           % Time step in the time grid
global X   % Vector [x0, x1, ..., xn] specifying the points at which a numerical solution is requested for every value in tspan ds
global T simT rawT   % Vector [t0, t1,..., tf] specifying the points at which a solution is requested for every value in X 
                     % T - what is fed to PDE solver
                     % simT, rawT - holders for simulation and experimental
                     % time grids respectively
                     
% Initialisations
simD1 = 10;   
simD2 = 0.5;  
simK1 = 0.1;  
simK2 = 0.01; 

fittedD1 = 10;   
fittedD2 = 0.5;  
fittedK1 = 0.1;  
fittedK2 = 0.01; 

roiSize = 10;
modelLength = 60;
deltaX = 1;
simulationTime = 20;
deltaT = 0.5;
noise1 = 0.002;
noise2 = 0.02;

% Loading GUI file 'fitgui.fig' 
fig_handle = fitgui;






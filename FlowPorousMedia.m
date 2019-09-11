function  FlowPorousMedia()
clc;
clear;
close all;
warning off;

% Import mesh and convert to second order mesh
[p,e,t] = importMeshGmsh('heatExchanger.msh');
[p,e,t,nVnodes] = convertMeshToSecondOrder(p,e,t);
%%
figure('Renderer','Opengl','Position',[10,10,720,640])
displayMesh2D(p,t);

%% Generate random permeability field of sand
% Define water dynamic viscosity [Pa*s]
mu = 0.00025;
K = (1e-9)*rand(nVnodes,1,'single');
figure('Position',[10,10,700,640])
displaySolution2D(p,t,K,'Porous media permeability');
colorbar off

%% Assemble global problem matrix accounting for Darcy law
[M, F] = assembleDiffusionMatrix2D(p,t,K);

% Set high water pressure head at inlet (rho*g*h) with h about 5.5 m
[M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,8,'value',500*9.81*7.5);

% Set zero pressure at outlet
[M,F] = imposeScalarBoundaryCondition2D(p,e,M,F,9,'value',0);

% Solve for pressure field
pressure = M\F;

% Display pressure field
figure('Position',[10,10,720,640])
displaySolution2D(p,t,pressure,'Pressure field');
colorbar off
%% Compute horizontal vx and vertical vy velocity
[px,py] = gradient(t);%solutionGradient2D(p,t,pressure);
vx = (1-K(4)/mu.*px);
vy = (1-K(6)/mu.*py);

% Compute velocity magnitude
vmag = sqrt((vx.^2+vy.^2)./(px.^2+py.^2)); 

% Display velocity field
figure('Renderer','Opengl','Position',[10,10,720,640])
displaySolution2D(p,t,vmag,'Velocity magnitude [m/s]');
exportToGmsh2D('out.msh', vmag, p, t,'ASCII');
colorbar off
drawnow;

end


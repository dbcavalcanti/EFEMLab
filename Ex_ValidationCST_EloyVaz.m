%% ========================= FINITE ELEMENT ===============================
%
% This script presents a simple example used to validate the implementation
% of the CST element. The example consists in a bar submitted to a traction
% modeled using four CST elements.
% 
% Reference: 
%   L. E. Vaz. Método dos Elementos Finitos em Análise de Estruturas.
%   Elsevier, 2011.
%
% Author: Danilo Cavalcanti
%
% Date: January, 27th, 2023.
%
%% ========================================================================
%
close all;
% Clear the classes to avoid using a non-updated version
clear Element RegularElement  EnrichedElement;
clear EnrichedElement_SOS EnrichedElement_KOS EnrichedElement_KSON;
clear Fracture Fracture_ConstantJump Fracture_LinearJump;
clear Shape Shape_CST Shape_ISOQ4;
clear Model IntPoint Result;
% Clear the workspace and the command window
clear; clc;
%Use all folders and subfolders
addpath(genpath('./')); 
%
%% ========================== MODEL CREATION ==============================
% --- Mesh of continuum elements ------------------------------------------

% Nodes' coordinates (mm)
NODE = [ 0.0   0.0;
         0.0   1.0;
         2.0   0.0;
         2.0   1.0
         4.0   0.0;
         4.0   1.0];

% Type of elements
type = 'CST';

% Element's connectivity
ELEM = [ 1 4 2;
         1 3 4;
         3 6 4;
         3 5 6];

% Thickness (mm)
t = 1.0;

% --- Mesh of the fracture elements ---------------------------------------

% Coordinates of the nodes that define the discontinuities (mm)
NODE_D = [];

% Fractures definition (by segments)
FRACT = [];

% --- Material properties of the domain -----------------------------------

% Define the material model of the continuum: 'elastic'
matModel = 'elastic';

% Material parameters
E  = 20000.0;          % Young's modulus (MPa)
nu = 0.2;             % Poisson's ratio
mat = [E  nu];    % Material parameters vector

% --- Material properties of the fracture ---------------------------------

% Define the traction constitutive law: 'elastic', 'isotropicDamage'
tractionLaw = 'elastic';  

kn = 0.0;           % Normal stiffness (MPa/mm)
ks = 0.0;           % Shear stiffness (MPa/mm)
matfract = [ks, kn];

% --- Analysis model ------------------------------------------------------

% Type of analysis: 'PlaneStress' or 'PlaneStrain'
anm = 'PlaneStress';

% --- Boundary conditions --------------------------------------------------

% Define supports
SUPP = zeros(size(NODE,1),2);
SUPP([1 2],:) = [1 1;1 0]; 

% Define prescribe displacements
PRESCDISPL = zeros(size(NODE,1),2);

% Define the load conditions
LOAD = zeros(size(NODE,1),2);
LOAD([5 6],:) = [10.0  0.0; 10.0  0.0];

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
intOrder = 1;

%% ========================= INITIALIZATION ===============================

% Create the model object
mdl = Model(NODE, ELEM, [], [], t, matModel, mat, tractionLaw, [], [], ...
            anm, type, SUPP, LOAD, PRESCDISPL, intOrder,'',...
            [], [],'', [], '', 0, zeros(size(ELEM,1)));

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Plot the mesh with the supports
mdl.plotMeshWithBC();

%% ========================== RUN ANALYSIS ================================

% Solve the structural analysis problem
anl = Anl_Linear();
% anl = Anl_Nonlinear();
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();
mdl.plotDeformedMesh();

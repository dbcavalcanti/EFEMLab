%% ==================== EMBEDDED FINITE ELEMENT ===========================
%
% Author: Danilo Cavalcanti
%
% Date: January, 27th, 2023.
%
%% ========================================================================
%
close all;
% Clear the classes to avoid using a non-updated version
clear Element EnrichedElement EnrichedElement_KOS EnrichedElement_KSON;
clear Fracture RegularElement Shape Shape_CST Shape_ISOQ4 Model Result;
% Clear the workspace and the command window
clear; clc;
%Use all folders and subfolders
addpath(genpath('./')); 
%
%% ========================== MODEL CREATION ==============================
% --- Mesh of continuum elements ------------------------------------------

% Nodes' coordinates (mm)
NODE = [  0.0     0.0;
        200.0     0.0;
        200.0   200.0;
          0.0   200.0];

% Type of elements
type = 'ISOQ4';

% Element's connectivity
ELEM = [1 2 3 4];

% Thickness (mm)
t = 1.0;

% --- Mesh of the fracture elements ---------------------------------------

% Coordinates of the nodes that define the discontinuities (mm)
NODE_D = [100.0    0.0;
          100.0  200.0];

% Fractures definition (by segments)
FRACT = [1 2];

% --- Material properties of the domain -----------------------------------

% Define the material model of the continuum: 'elastic'
matModel = 'elastic';

% Material parameters
E   = 30.0e3;     % Young's modulus (MPa)
nu  = 0.0;        % Poisson's ratio
mat = [E  nu];    % Material parameters vector

% --- Material properties of the fracture ---------------------------------

% Define the traction constitutive law: 'elastic', 'isotropicDamage'
tractionLaw = 'elastic';  

% Flag to apply a penalization on compression 
tractionLawPenal = true;

% Values of the material constitutive model parameters
kn = 0.0;           % Normal stiffness (MPa/mm)
ks = 0.0;           % Shear stiffness (MPa/mm)
matfract = [ks, kn];

% --- Analysis model ------------------------------------------------------

% Type of analysis: 'PlaneStress' or 'PlaneStrain'
anm = 'PlaneStress';

% --- Boundary conditions --------------------------------------------------

% Define supports
SUPP = zeros(size(NODE,1),2);
SUPP([1 2 3 4],:) = [1 1;1 1;1 1;1 1];

% Define prescribe displacements
PRESCDISPL = zeros(size(NODE,1),2);
displ = 10.0;
PRESCDISPL([2 3],:) = [0.0   -displ;
                       0.0    displ];

% Define the load conditions
LOAD = zeros(size(NODE,1),2);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
intOrder = 2;

%% ===================== EFEM FORMULATION SETUP ===========================

% Type of formulation
enhancementType = 'KSON';

% Apply a sub-division of the domain to perform the numerical integration
subDivInt = true;

% Consider the stretch part of the mapping matrix
stretch = true;

% Order of the interpolation of the jump displacement field
jumpOrder = 1;

%% ========================= PRE-PROCESSING ===============================

% Compute the matrix to identify to which element the fracture belongs
IDenr = 1; 

%% ========================= INITIALIZATION ===============================

% Create the model object
mdl = Model(NODE, ELEM, NODE_D, FRACT, t, matModel, mat, tractionLaw, ...
            tractionLawPenal, matfract, anm, type, SUPP, LOAD, ...
            PRESCDISPL, intOrder, enhancementType, subDivInt, stretch, ...
            jumpOrder, IDenr);

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Plot the mesh with the supports
mdl.plotMeshWithBC();

%% ========================== RUN ANALYSIS ================================

% Solve the structural analysis problem
anl = Anl_Linear();
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();
mdl.plotDeformedMesh();
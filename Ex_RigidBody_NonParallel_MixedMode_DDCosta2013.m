%% ==================== EMBEDDED FINITE ELEMENT ===========================
%
% This script consists in the example presented in the paper by
% Dias-da-Costa et al. (2009), in the section 4.5.1.
% This exampel simulates a rigid body motion crack opening. The crack is
% parallel to the element border and a mode I is induced.
%
%
% Reference:
% DIAS DA COSTA, D.; ALFAIATE, J.; SLUYS, L. J.; JÚLIO, E. Towards a
% generalization of a discrete strong discontinuity approach. Computer 
% Methods in Applied Mechanics and Engineering, v. 198, n. 47-48, 
% p. 3670–3681, oct 2009.
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
NODE = [0.0   0.0;
        2.0   0.0;
        2.0   2.0;
        0.0   2.0];

% Type of elements
type = 'ISOQ4';

% Element's connectivity
ELEM = [1 2 3 4];

% Thickness (mm)
t = 1.0;

% --- Mesh of the fracture elements ---------------------------------------

% Coordinates of the nodes that define the discontinuities (mm)
NODE_D = [0.0  0.25;
          2.0  1.75];

% Fractures definition (by segments)
FRACT = [1 2];

% --- Material properties of the domain -----------------------------------

% Define the material model of the continuum: 'elastic'
matModel = 'elastic';

% Material parameters
E   = 1.0e8;      % Young's modulus (MPa)
nu  = 0.0;        % Poisson's ratio
mat = [E  nu];    % Material parameters vector

% --- Material properties of the fracture ---------------------------------

% Define the traction constitutive law: 'elastic', 'isotropicDamage'
tractionLaw = 'elastic';  

% Flag to apply a penalization on compression 
tractionLawPenal = true;

% Values of the material constitutive model parameters
kn   = 1.0e0;            % Normal stiffness (MPa/mm)
ks   = 1.0e0;            % Shear stiffness (MPa/mm)
ft   = 3.0;              % Tensile strength
Gf   = 1.0;              % Fracture energy
beta = 0.0;              % Shear factor

% Assemble the vector with the material properties
matfract = [ks, kn, ft, Gf];

% --- Analysis model ------------------------------------------------------

% Type of analysis: 'PlaneStress' or 'PlaneStrain'
anm = 'PlaneStress';

% --- Boundary conditions --------------------------------------------------

% Define supports
SUPP = zeros(size(NODE,1),2);
SUPP([1 2],:) = [1 1;1 1];

% Define prescribe displacements
PRESCDISPL = zeros(size(NODE,1),2);

% Define the load conditions
LOAD = zeros(size(NODE,1),2);
LOAD(4,:) = [-0.5 1.5]; 

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
intOrder = 2;

%% ===================== EFEM FORMULATION SETUP ===========================

% Type of formulation
enhancementType = 'KOS';

% Apply a sub-division of the domain to perform the numerical integration
subDivInt = false;

% Consider the stretch part of the mapping matrix
stretch = false;

% Order of the interpolation of the jump displacement field
jumpOrder = 1;

% Enrichment degree of freedom ('w' or 'alpha')
enrVar = 'w';

%% ========================= PRE-PROCESSING ===============================

% Compute the matrix to identify to which element the fracture belongs
IDenr = 1; 

%% ========================= INITIALIZATION ===============================

% Create the model object
mdl = Model(NODE, ELEM, NODE_D, FRACT, t, matModel, mat, tractionLaw, ...
            tractionLawPenal, matfract, anm, type, SUPP, LOAD, ...
            PRESCDISPL, intOrder, enhancementType, subDivInt, stretch, ...
            enrVar, jumpOrder, IDenr);

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

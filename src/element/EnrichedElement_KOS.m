%% EnrichedElement_KOS class
%
% This class implements a finite element with an embedded strong 
% discontinuity using the KOS formulation.
%
% Reference of the element formulation:
%
% Dias-da-Costa, D. Alfaiate, J. Sluys, L. J. Júlio, E. A discrete strong
% discontinuity approach. Engineering Mechanics, vol. 76. 2009.
%
% Reference of the classification:
%
% Jirasek, M. Comparative study on finite elements with embedded
% discontinuities. Computat. Methods Appl. Mech. Engrg., vol 188. 2000.
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: January, 2023
%
%% Class definition
classdef EnrichedElement_KOS < EnrichedElement
    %% Constructor method
    methods
        function this = EnrichedElement_KOS(type, node, elem, anm, t, matModel, mat, nGP, gla, fracture, glw, subDivInt, stretch, enrVar, jumpOrder, staticCondensation)
            this = this@EnrichedElement(type, node, elem, anm, t, matModel, mat, nGP, gla, fracture, glw, subDivInt, stretch, enrVar, jumpOrder, staticCondensation);
        end
    end
    %% Public methods
    % Implementation of the abstract methods declared in super-class
    methods

        % -----------------------------------------------------------------
        % Compute the matrix Gr. This matrix discretized the bounded part
        % of the enhanced strains wrt to the vector of enhanced dof (w).
        function [Gr, Gv] = enhancedStrainCompatibilityMtrcs(this, B, Xn)

            % Compute the matrix Gr
            Gr = this.enhancedKinematicCompatibilityMtrx(B, Xn);

            % Compute the matrix Gv
            Gv = Gr;

        end
        
    end

end
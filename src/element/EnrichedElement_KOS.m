%% EnrichedElement_KOS class
%
% This class implements an ISO-Q4 with an embedded strong discontinuity
% using the KOS formulation.
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
        % This matrix is being defined considering the simplification that
        % only rigid-body displacements are incorporated in the jum
        % displacement field.
        function Gc = enhancedStrainCompatibilityMtrx(this, B, Xn)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);

            % Evaluate the Heaviside function at the point X
            h = this.heavisideFnc(X);

            % Compute the Heaviside matrix
            Hd = this.heavisideMtrx();

            % Compute the mapping matrix
            M = this.elementMappingMtrx();

            % Identity matrix
            I = eye(size(Hd,1));

            % Compute the matrix Gr
            Gc = B*(h*I - Hd)*M;

        end
            
        % -----------------------------------------------------------------
        % Compute the matrix Gv. This matrix discretized the bounded part
        % of the virtual enhanced strains wrt to the vector of enhanced
        % dof (w). The KOS formulation (see Jirasek (2000) is considered),
        % where Gv = Gr.
        function Gv = enhancedStressCompatibilityMtrx(this, B, Xn)

             Gv = this.enhancedStrainCompatibilityMtrx(B, Xn);

        end

    end

end
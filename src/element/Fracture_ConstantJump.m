%% Fracture_ConstantJump class
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: February, 2023
%
%% Class definition
classdef Fracture_ConstantJump < Fracture
    %% Constructor method
    methods
        function this = Fracture_ConstantJump(node, elem, t, matModel, mat, glw, penal)
            this = this@Fracture(node, elem, t, matModel, mat, glw, penal);
        end
    end
    %% Public methods
    % Implementation of the abstract methods declared in super-class
    methods

        %------------------------------------------------------------------
        % This function computes the matrix of the shape function
        % to evaluate the displacement jump based on the enrichment degrees
        % of freedom 'alpha'.
        function N = interpJumpShapeMtrx(~,~,~)

            % Shape function matrix
            N = [ 1.0  0.0 ;
                  0.0  1.0 ];

        end

        % -----------------------------------------------------------------
        % Compute the jump transmission matrix M. This matrix relates the
        % enrichment degrees of freedom alpha with the enhanced
        % displacement field.
        function M = jumpTransmissionMtrx(~,~,~,~,~)

            M = [ 1.0  0.0;
                  0.0  1.0 ];

        end

        %------------------------------------------------------------------
        % This function computes the element's rotation matrix. Change from
        % the local coordinate system mn to the global system xy
        function R = rotationMtrx(this)

            % Rotation of a point
            R = [ this.m(1)   this.m(2);
                  this.n(1)   this.n(2) ];

        end

        % -----------------------------------------------------------------
        % Matrix to transform the enrichment degrees of freedom from alpha
        % to w
        function Se = transformAlphaToW(~)

            % Matrix Se
            Se = [ 1.0  0.0;
                   0.0  1.0];

        end

        % -----------------------------------------------------------------
        % Matrix to transform the enrichment degrees of freedom from w
        % to alpha
        function S = transformWToAlpha(~)

            % Matrix S
            S = [ 1.0  0.0;
                  0.0  1.0];

        end 

        %------------------------------------------------------------------
        % This function compute the stress interpolation vector
        function S = stressIntVct(this, shape, node)
            S = this.stressIntVctFnc(shape, node, 0);
        end
 
    end

end
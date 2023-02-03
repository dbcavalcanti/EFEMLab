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
classdef Fracture_LinearJump < Fracture
    %% Constructor method
    methods
        function this = Fracture_LinearJump(node, elem, t, matModel, mat, glw, penal)
            this = this@Fracture(node, elem, t, matModel, mat, glw, penal);
        end
    end
    %% Public methods
    % Implementation of the abstract methods declared in super-class
    methods

        % -----------------------------------------------------------------
        % Compute the jump transmission matrix M. This matrix relates the
        % enrichment degrees of freedom alpha with the enhanced
        % displacement field.
        function M = jumpTransmissionMtrx(this,X,enrVar,stretch,nu)

            % Initialize the mapping matrix
            M = zeros(2,4);

            % Add the translation transmission
            M = this.addTranslationMtrx(M);

            % Add the relative rigid-body rotation transmission
            M = this.addRelativeRotationMtrx(X,M);

            % Add the stretching along the tangential direction
            % transmission
            if stretch(1) == true
                M = this.addStretchingMtrx(X,M);
    
                % Add the stretching along the normal direction due to the
                % Poisson effect transmission
                if stretch(2) == true
                    M = this.addPoissonEffectMtrx(X,M,nu);
                end
            end

            % Transform the enrichment variables from alpha to w
            if strcmp(enrVar,'w')
                Se = this.transformAlphaToW();
                M  = M*Se;
            end

        end
 
    end

    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Mapping matrix of the translation part of the jump displacement
        % field 
        function M = addTranslationMtrx(~,M)

            % Translation matrix
            Mt = [1.0   0.0   0.0   0.0;
                  0.0   1.0   0.0   0.0]; 

            % Add to the jump matrix
            M = M + Mt;

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the relative rotation part of the jump 
        % displacement field. Evaluated at one point X.
        function M = addRelativeRotationMtrx(this,X,M)
            
            % Fracture geometric properties
            m    = this.m;
            cs   = m(1);
            sn   = m(2);
            Xref = this.Xref;
            
            % Relative position vector
            DX = X - Xref;

            % Relative rotation transmission matrix matrix
            Mrot = [ 0.0  0.0   sn*DX(2)  -cs*DX(2);
                     0.0  0.0  -sn*DX(1)   cs*DX(1)];

            % Add to the jump matrix
            M = M + Mrot;

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the stretching part of the jump displacement 
        % field. Evaluated at one point X.
        function M = addStretchingMtrx(this,X,M)
            
            % Fracture geometric properties
            m    = this.m;
            cs   = m(1);
            sn   = m(2);
            Xref = this.Xref;

            % Tangential relative coordinate
            s = m*(X' - Xref');

            % Mapping matrix
            Ms = [ 0.0  0.0  s*cs*cs  s*cs*sn;
                   0.0  0.0  s*cs*sn  s*sn*sn];

            % Add to the jump matrix
            M = M + Ms;

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the stretching part in the normal direction due
        % to the Poisson effect. Evaluated at one point X.
        function M = addPoissonEffectMtrx(this,X,M,nu)
            
            % Fracture geometric properties
            m    = this.m;
            cs   = m(1);
            sn   = m(2);
            Xref = this.Xref;

            % Relative position vector
            DX = X - Xref;

            % Auxiliary coefficients
            a1 = sn*sn*DX(1) - cs*sn*DX(2);
            a2 = cs*cs*DX(2) - cs*sn*DX(1);

            % Mapping matrix
            Mnu = [ 0.0  0.0  -nu*cs*a1  -nu*sn*a1;
                    0.0  0.0  -nu*cs*a2  -nu*sn*a2];

            % Add to the jump matrix
            M = M + Mnu;

        end

        % -----------------------------------------------------------------
        % Matrix to transform the enrichment degrees of freedom from alpha
        % to w
        function Se = transformAlphaToW(this)

            % Find the tangential coordinates along the discontinuity of
            % the initial and final nodes
            s1 = this.tangentialLocCoordinate(0.0);
            s2 = this.tangentialLocCoordinate(1.0);

            % Matrix Se
            Se = [-s2    0.0   s1    0.0;
                   0.0   -s2   0.0    s1;
                   1.0   0.0  -1.0   0.0;
                   0.0   1.0   0.0  -1.0]/(s1 - s2);

        end

        % -----------------------------------------------------------------
        % Matrix to transform the enrichment degrees of freedom from w
        % to alpha
        function S = transformWToAlpha(this)

            % Find the tangential coordinates along the discontinuity of
            % the initial and final nodes
            s1 = this.tangentialLocCoordinate(0.0);
            s2 = this.tangentialLocCoordinate(1.0);

            % Matrix Se
            S = [ 1.0   0.0    s1   0.0;
                  0.0   1.0   0.0    s1;
                  1.0   0.0    s2   0.0;
                  0.0   1.0   0.0    s2]/(s1 - s2);

        end 


    end

end
%% EnrichedElement_KSON class
%
% This class implements an ISO-Q4 with an embedded strong discontinuity
% using the KSON formulation.
%
% Reference of the element formulation:
%
% Linder, C. Armero, F. Finite elements with embedded strong
% discontinuities for the modeling of failure in solids. 
% Int. J. Numer. Meth. Engng., vol. 72. 2007.
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
classdef EnrichedElement_KSON < EnrichedElement
    %% Constructor method
    methods
        function this = EnrichedElement_KSON(type, node, elem, t, mat, nGP, gla, fracture, glw, subDivInt, stretch, jumpOrder)
            this = this@EnrichedElement(type, node, elem, t, mat, nGP, gla, fracture, glw, subDivInt, stretch, jumpOrder);
        end
    end
    %% Public methods
    % Implementation of the abstract methods declared in super-class
    methods

        % -----------------------------------------------------------------
        % Compute the matrix Gc. This matrix discretized the bounded part
        % of the enhanced strains wrt to the vector of enhanced dof (w).
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

            % Compute the matrix Gc
            Gc = B*(h*I - Hd)*M;

        end
            
        % -----------------------------------------------------------------
        % Compute the matrix Gv. This matrix discretized the bounded part
        % of the virtual enhanced strains wrt to the vector of enhanced
        % dof (w). 
        function Gv = enhancedStressCompatibilityMtrx(this, ~, Xn)

            % Point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node,Xn);

            % Centroid of the element
            X0 = this.shape.coordNaturalToCartesian(this.node,[0,0]);

            % Get the fracture normal vector
            n = this.fracture.n';

            % Get the fracture length
            ld = this.fracture.ld;

            % Compute the projection matrix
            P = [n(1) 0.0;
                 0.0  n(2);
                 n(2) n(1)];
            R = [n(1) n(2);
                -n(2) n(1)];
            P = P*R';

            % Compute the coefficients of the polynomial interpolation
            C0 = getPolynomialCoeffs(this,X0,ld,0);
            C1 = getPolynomialCoeffs(this,X0,ld,1);

            % Compute the cartesian coordinates in the element's local
            % system
            Xrel = X - X0;

            % Evaluate the polynomial
            g0 = C0(1) + C0(2)*Xrel(1) + C0(3)*Xrel(2);
            g1 = C1(1) + C1(2)*Xrel(1) + C1(3)*Xrel(2);

            % Compute the columns of Gv:
            G0 = -g0*P;
            G1 = -g1*P;

            % Assemble the operator Gv
            Gv = [G0, G1];

            % Transform for considering the nodal jumps vector [w] as
            % enhancement dofs
            Se = this.transformAlphaToW();
            Gv = Gv*Se;

        end

        % -----------------------------------------------------------------
        % Compute the coefficients of the polynomial g^(k) that is used to 
        % approximate the submatrix matrix Gv^(k).
        % The submatrix is defined as:
        %   Gv^(k) = g^(k) * P. 
        % Where:
        %   g^(k) = c0 + c1 * xrel + c2 * yrel
        % 
        function C = getPolynomialCoeffs(this,X0,ld,k)

            % Get integration points
            [Xip,w,nIntPoints] = this.shape.getIntegrationPoints(this.intOrder);
            
            % Numerical integration of the stiffness matrix components
            Gramm = zeros(3,3);
            for i = 1:nIntPoints

                % Integration point in the global cartesian coordinate
                % system
                [~, detJ] = this.shape.BMatrix(this.node,[Xip(1,i), Xip(2,i)]);
                X = this.shape.coordNaturalToCartesian(this.node,[Xip(1,i), Xip(2,i)]);

                % Integration point in the element local cartesian
                % coordinate system, with origin located at its
                % centroid
                Xrel = X - X0;
        
                % Compute the elastic constitutive matrix
                dA = [  1.0          Xrel(1)           Xrel(2);
                       Xrel(1)    Xrel(1)*Xrel(1)   Xrel(2)*Xrel(1);
                       Xrel(2)    Xrel(1)*Xrel(2)   Xrel(2)*Xrel(2)];
        
                % Numerical integration coefficient
                c = w(i)*detJ;
        
                % Numerical integration of the stiffness matrix Kaa
                Gramm = Gramm + dA * c;
            
            end

            % Compute the integration point of the discontinuity
            [w,gp] = this.fracture.getlineQuadrature_GaussLobatto();

            % Numerical integration
            vec = zeros(3,1);
            for i = 1:length(gp)
                WdetJ = w(i)*ld;
                Nw    = this.fracture.shapeFncFracture(gp(i));
                X     = Nw*[this.fracture.node(1,:),this.fracture.node(2,:)]';
                Xrel = X - X0;
                s     = this.fracture.m*(X - this.fracture.Xref');
                dVec  = [s^k; (s^k)*Xrel(1); (s^k)*Xrel(2)];
                vec   = vec + dVec * WdetJ;
            end

            % Compute the coefficients
            C = Gramm \ vec;

        end

        function Se = transformAlphaToW(this)
            Nw    = this.fracture.shapeFncFracture(0);
            X     = Nw*[this.fracture.node(1,:),this.fracture.node(2,:)]';
            s1     = this.fracture.m*(X - this.fracture.Xref');

            Nw    = this.fracture.shapeFncFracture(1);
            X     = Nw*[this.fracture.node(1,:),this.fracture.node(2,:)]';
            s2     = this.fracture.m*(X - this.fracture.Xref');

            Se = [-s2    0   s1    0;
                    0  -s2    0   s1;
                    1    0   -1    0;
                    0    1    0   -1]/(s1 - s2);

        end

    end

end
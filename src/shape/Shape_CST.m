%% Shape_CST class
%
% This is class defines the behavior of a constant strain triangular
% isoparametric element.
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version:January 2023
%
classdef Shape_CST < Shape
    %% Constructor method
    methods
        function this = Shape_CST()
            this = this@Shape('CST');
        end
    end

    %% Public methods: methods defined in the abstract superclass
    methods

        % -----------------------------------------------------------------
        % Evaluate the shape function at a given point X of a linear 
        % quadrilateral isoparametric element.
         function N = shapeFnc(~,Xn)

            % Natural coordinates of the given point
            r = Xn(1); s = Xn(2);

            % Shape functions
            N1 = 1.0 - r - s;
            N2 = r;
            N3 = s;
            N   = [ N1  N2  N3 ];

         end

         % -----------------------------------------------------------------
         % Get the shape function matrix
         function Nm = shapeFncMtrx(this,Xn)

             % Vector with the shape functions
             N = this.shapeFnc(Xn);
            
             % Shape function matrix
             Nm = [N(1)  0.0   N(2)  0.0   N(3)  0.0;
                   0.0   N(1)  0.0   N(2)  0.0   N(3)];

         end

         % -----------------------------------------------------------------
         % Compute the derivatives of the shape functions wrt to the
         % natural coordinate s
         function dNdxn = shapeFncDrv(~,~)

            % Derivatives of the shape functions
            dN1_dr = -1.0;    dN1_ds = -1.0;
            dN2_dr =  1.0;    dN2_ds =  0.0;
            dN3_dr =  0.0;    dN3_ds =  1.0;
            dNdxn   = [ dN1_dr  dN2_dr  dN3_dr;
                        dN1_ds  dN2_ds  dN3_ds];

         end

         % -----------------------------------------------------------------
         % Compute the jacobian matrix
         function J = JacobianMtrx(this,X,Xn)

            % Compute the shape function derivatives wrt to the natural
            % coordinate system
            dNdxn = this.shapeFncDrv(Xn);
            
            % Jacobian matrix
            J = dNdxn*X;

         end

         % -----------------------------------------------------------------
         % Compute the determinant of the jacobian
         function detJ = detJacobian(this,X,Xn)
              
            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);
            detJ = det(J);

         end

         % -----------------------------------------------------------------
         % Compute the strain-displacement matrix
         function [B,detJ] = BMatrix(this,X,Xn)

            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);

            % Determinant of the Jacobian matrix
            detJ = det(J);

            % Compute the derivatives of the shape functions wrt to the
            % natural coordinate system
            dNdxn = this.shapeFncDrv(Xn);

            % Compute the derivatives of the shape functions wrt to the
            % global cartesian coordinate system
            dNdx = J\dNdxn;

            B = zeros(3,3*2);
            for i = 1:3
                B(1,2*i-1) = dNdx(1,i); 
                B(2,2*i)   = dNdx(2,i);
                B(3,2*i-1) = dNdx(2,i); B(3,2*i) = dNdx(1,i);
            end

         end

         % -----------------------------------------------------------------
         % Get the integration points:
         % Output:
         %      X: Coordinates of the integration points in the natural
         %         system
         %      w: Weight associated to each integration point
         %      nIntPoints: Total number of integration points
         function [X,w,nIntPoints] = getIntegrationPoints(~,order,~)

             if order == 2
                 X          = [1/6,  1/6;
                               2/3,  1/6;
                               1/6,  2/3]';  
                 w          = [1/6,1/6,1/6];
                 nIntPoints = 3;
             elseif order == 1
                 X          = [1/3, 1/3]';
                 w          = 0.5;
                 nIntPoints = 1;
             end

         end

         % -----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of an
         %           CST element
         %   Xn   : vector with the xi and eta coordinates of a point in 
         %          the natural coordinate system
         %
         % Output:
         %   X : vector with the x and y coordinates of a point in the 
         %       global coordinate system
         function X = coordNaturalToCartesian(this,NODE,Xn)

            % Extract the nodal coordinates
            x = NODE(:,1);
            y = NODE(:,2);
            
            % Vector with the shape functions
            Nv = this.shapeFnc(Xn);
            
            % Initialize output
            X = [0.0, 0.0];
            
            % Interpolation the position
            X(1) = Nv(1)*x(1) +  Nv(2)*x(2) +  Nv(3)*x(3);
            X(2) = Nv(1)*y(1) +  Nv(2)*y(2) +  Nv(3)*y(3);

         end

         % ----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Reference: Section 7.3 from the book Concepts and Applications
         % of Finite Element Analysis from Cook et al. (2001)
         function Xn = coordCartesianToNatural(this,NODE,X)
            
            % Compute the area of the triangle
            A = this.areaTriangle(NODE(1,:), NODE(2,:), NODE(3,:));

            % Compute the areas of the sub-triangles
            A1 = this.areaTriangle(NODE(2,:), NODE(3,:), X);
            A2 = this.areaTriangle(NODE(3,:), NODE(1,:), X);

            % Compute the area coordinates
            xi1 = A1/A;
            xi2 = A2/A;
            xi3 = 1.0 - xi1 - xi2;

            Xn = [xi2, xi3];
            
         end

        % -----------------------------------------------------------------
        % Compute the area of a triangle given the three point in
        % counterclock-wise order
        function A = areaTriangle(~, P1, P2, P3)
            P1P2 = P2 - P1;
            P1P3 = P3 - P1;
            A = (P1P2(1)*P1P3(2) - P1P2(2)*P1P3(1))*0.5;
        end

    end

    methods(Static)
        function [X,w,nIntPoints] = getTriangleLinearQuadrature()
            % Coordinates of the integration points
            X          = [1/6,  1/6;
                          2/3,  1/6;
                          1/6,  2/3]';  
             % Weights
             w          = [1/6,1/6,1/6];
             % Total number of integration points
             nIntPoints = 3;
        end

        
    end
end
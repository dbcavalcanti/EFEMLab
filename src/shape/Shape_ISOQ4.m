%% Shape_ISOQ4 class
%
% This is class defines the behavior of a linear quadrilateral
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
classdef Shape_ISOQ4 < Shape
    %% Constructor method
    methods
        function this = Shape_ISOQ4()
            this = this@Shape('ISOQ4');
        end
    end

    %% Public methods: methods defined in the abstract superclass
    methods

        % -----------------------------------------------------------------
        % Evaluate the shape function at a given point X of a linear 
        % quadrilateral isoparametric element.
         function N = shapeFnc(~,Xn)

            % Natural coordinates of the given point
            xi = Xn(1); eta = Xn(2);

            % Shape functions
            N1 = (1.0 - xi)*(1.0 - eta)/4.0;
            N2 = (1.0 + xi)*(1.0 - eta)/4.0;
            N3 = (1.0 + xi)*(1.0 + eta)/4.0;
            N4 = (1.0 - xi)*(1.0 + eta)/4.0;
            N   = [ N1  N2  N3  N4 ];

         end

         % -----------------------------------------------------------------
         % Get the shape function matrix
         function Nm = shapeFncMtrx(this,Xn)

             % Vector with the shape functions
             N = this.shapeFnc(Xn);
            
             % Shape function matrix
             Nm = [N(1)  0.0   N(2)  0.0   N(3)  0.0   N(4)  0.0;
                   0.0   N(1)  0.0   N(2)  0.0   N(3)  0.0   N(4)];

         end

         % -----------------------------------------------------------------
         % Compute the derivatives of the shape functions wrt to the
         % natural coordinate s
         function dNdxn = shapeFncDrv(~,Xn)

            % Natural coordinates of the given point
            xi = Xn(1); eta = Xn(2);

            % Derivatives of the shape functions
            dN1_dxi = -(1-eta)/4;    dN1_deta = -(1-xi)/4;
            dN2_dxi =  (1-eta)/4;    dN2_deta = -(1+xi)/4;
            dN3_dxi =  (1+eta)/4;    dN3_deta =  (1+xi)/4;
            dN4_dxi = -(1+eta)/4;    dN4_deta =  (1-xi)/4;
            dNdxn   = [ dN1_dxi   dN2_dxi   dN3_dxi   dN4_dxi ;
                        dN1_deta  dN2_deta  dN3_deta  dN4_deta ];

         end

         % -----------------------------------------------------------------
         % Compute the jacobian matrix
         function J = JacobianMtrx(this,X,Xn)

            % Compute the shape function derivatives wrt to the natural
            % coordinate system
            dNdxn = this.shapeFncDrv(Xn);
              
            % Jacobian matrix
            J = dNdxn * X;

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

            B = zeros(3,4*2);
            for i = 1:4
                B(1,2*i-1) = dNdx(1,i); 
                B(2,2*i)   = dNdx(2,i);
                B(3,2*i-1) = dNdx(2,i); B(3,2*i) = dNdx(1,i);
            end

         end

         % -----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of an
         %           ISOQ4 element
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
            X(1) = Nv(1)*x(1) +  Nv(2)*x(2) +  Nv(3)*x(3) +  Nv(4)*x(4);
            X(2) = Nv(1)*y(1) +  Nv(2)*y(2) +  Nv(3)*y(3) +  Nv(4)*y(4);

         end

         % -----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of an 
         %          ISOQ4 element
         %   X    : vector with the x and y coordinates of a point in the 
         %          global cartesian coordinate system
         %
         % Output:
         %   xi, eta: coordinates xi and eta of the point X in the natural
         %            coordinate system.
         %
         % Reference:
         %   Felippa, C. A. Introduction to the Finite Element Methods. 
         %   2004.(Chapter 23, item 6)
         function Xn = coordCartesianToNatural(~,NODE,X)
            
            
            % Extract the nodal coordinates
            x = NODE(:,1);
            y = NODE(:,2);
            
            % Felippa's algorithm:
            xb  = x(1) - x(2) + x(3) - x(4);
            yb  = y(1) - y(2) + y(3) - y(4);
            xcx = x(1) + x(2) - x(3) - x(4);
            ycx = y(1) + y(2) - y(3) - y(4);
            xce = x(1) - x(2) - x(3) + x(4);
            yce = y(1) - y(2) - y(3) + y(4);
            A   = ((x(3) - x(1))*(y(4) - y(2)) - (x(4) - x(2))*(y(3) - y(1))) / 2.0;
            J1  = (x(3) - x(4))*(y(1) - y(2)) - (x(1) - x(2))*(y(3) - y(4));
            J2  = (x(2) - x(3))*(y(1) - y(4)) - (x(1) - x(4))*(y(2) - y(3));
            x0  = (x(1) + x(2) + x(3) + x(4)) / 4.0;
            y0  = (y(1) + y(2) + y(3) + y(4)) / 4.0;
            xp0 = X(1) - x0;
            yp0 = X(2) - y0;
            bxi =  A - xp0*yb + yp0*xb;
            bet = -A - xp0*yb + yp0*xb;
            cxi = xp0*ycx - yp0*xcx;
            cet = xp0*yce - yp0*xce;
            
            % Natural coordinates
            xi  = 2.0*cxi / (-sqrt(bxi^2 - 2*J1*cxi) - bxi);
            eta = 2.0*cet / (sqrt(bet^2 + 2*J2*cet) - bet);
            Xn = [xi, eta];
            
         end

        %------------------------------------------------------------------
        % Function to get the matrices with the coordinates and weights of
        % integration points in the natural coordinate system. 
        % The definition of the integration points depends on the choice of
        % the order of the quadrature rule and the option for subdividing
        % the domain or not.
        % Output:
        %
        %   X : matrix with the coordinates of the integration points in
        %       the natural coordinate system. Each column of this matrix
        %       gives the coordinates of a point.
        %   W : vector with the weights associated with each point
        %   n : number of integration points
        %
        function [X,W,n] = getIntegrationPoints(this,intOrder,elem)
            if nargin == 2
                subDivInt = false;
            else
                subDivInt = elem.subDivInt;
            end

            if subDivInt == false
                % Compute the integration points
                [w,x] = this.getlineQuadrature(intOrder);

                % Assemble the matrix with the integration points 
                [Xip,Yip] = meshgrid(x,x);
                X = [Xip(:)';Yip(:)'];

                % Assemble the vector with the weights
                [Wxip,Wyip] = meshgrid(w,w);
                W = Wxip(:)' .* Wyip(:)';

                % Number of integration points
                n = size(X,2);

            else
                % Heaviside function evaluated in the nodes of the element
                h = elem.heavisideFnc(elem.node);

                % Get the nodes located at each region of the discontinuity
                Node_minus = [elem.node((h < 1),:); elem.fracture.node];
                Node_plus  = [elem.node((h > 0),:); elem.fracture.node];
                
                % Add the center of each region to the list of nodes 
                Node_minus = [Node_minus; mean(Node_minus(:,1))  mean(Node_minus(:,2))];
                Node_plus  = [Node_plus;  mean(Node_plus(:,1))   mean(Node_plus(:,2))];

                % Apply Delaunay's triangularization of each sub-region
                tri_minus = delaunay(Node_minus);
                tri_plus  = delaunay(Node_plus);

                % Get the integration points of a triangle
                shapeCST = Shape_CST();
                [x,w,nIntPoints] = shapeCST.getIntegrationPoints(intOrder,[]);

                % Total number of integration points
                n = nIntPoints*(size(tri_plus,1) + size(tri_minus,1));

                % Initialize the matrix of coordinates and the vector of
                % weights
                W = zeros(1,n);
                X = zeros(2,n);

                % Fill the matrix X and the vector W
                count = 1;
                for subelem = 1:size(tri_minus,1)
                    for ip = 1:length(w)
                        % Find the global cartesian coordinates of the integration point
                        [Xi] = shapeCST.coordNaturalToCartesian(...
                            Node_minus(tri_minus(subelem,:),:),x(:,ip));

                        % Weight
                        detJ = 1;%shapeCST.detJacobian(Node_minus(tri_minus(subelem,:),:),x(:,ip));
                        W(count) = w(ip)*detJ;
                
                        % Mapping the integration point to the quadrilateral
                        % element natural system
                        X(:,count) = this.coordCartesianToNatural(elem.node,Xi)';

                        count = count + 1;

                    end
                end
                for subelem = 1:size(tri_plus,1)
                    for ip = 1:length(w)
                        % Find the global cartesian coordinates of the integration point
                        [Xi] = shapeCST.coordNaturalToCartesian(...
                            Node_plus(tri_plus(subelem,:),:),x(:,ip));
                        
                        % Weight
                        detJ = 1;%shapeCST.detJacobian( Node_plus(tri_plus(subelem,:),:),x(:,ip));
                        W(count) = w(ip)*detJ;
                
                        % Mapping the integration point to the quadrilateral
                        % element natural system
                        X(:,count) = this.coordCartesianToNatural(elem.node,Xi)';

                        count = count + 1;

                    end
                end

            end

        end

        % -----------------------------------------------------------------
        % Integrand to compute the Gram Matrix
        function dH = integrandGramMtrx(this, node, X)

            X    = this.coordNaturalToCartesian(node,X);
            X0   = this.coordNaturalToCartesian(node,[0.0;0.0]);
            Xrel = X - X0;

            dH = [  1.0          Xrel(1)           Xrel(2);
                   Xrel(1)    Xrel(1)*Xrel(1)   Xrel(2)*Xrel(1);
                   Xrel(2)    Xrel(1)*Xrel(2)   Xrel(2)*Xrel(2)];
        end

        % -----------------------------------------------------------------
        % Integrand to compute the stress interpolation vector
        function n = getSizeStressIntVct(~)
            n = 3;
        end

        % -----------------------------------------------------------------
        % Integrand to compute the stress interpolation vector
        function dS = integrandStressIntVct(~,s,X,jumpOrder)

            X0   = this.coordNaturalToCartesian(node,[0.0;0.0]);
            Xrel = X - X0;

            if jumpOrder == 0
                dS = [  1.0;
                      Xrel(1);
                      Xrel(2)];
            elseif jumpOrder == 1
                dS = [  1.0       s;
                      Xrel(1)  s*Xrel(1);
                      Xrel(2)  s*Xrel(2)];
            end
        end

    end

    methods (Static)
        %------------------------------------------------------------------
        % Compute the Gauss integration points for a given order
        function [w,gp] = getlineQuadrature(order)
            % Number of Gauss points
            intOrder = order;
        
            % Define the coefficients of the three-term recurrence
            % relationship
            a = @(intOrder)(2*intOrder+1)./(intOrder+1);
            b = @(intOrder)0;
            c = @(intOrder)intOrder./(intOrder+1);
        
            % Constructe the symmetric tridiagonal matrix
            A = -b(0:intOrder-1)./a(0:intOrder-1); B = sqrt(c(1:intOrder-1)./(a(0:intOrder-2).*a(1:intOrder-1)));
            J = diag(B,1) + diag(A) + diag(B,-1);  
        
            % Compute the eigenvalues and eigenvectors
            [V,D] = eig(J,'vector');
        
            % Save (sorted) points and weights
            [gp,I] = sort(D);
            w = (2*V(1,I).^2)'; 
        
            % Note: The next three lines insure zero is zero and the points
            % and weights are perfectly symmetric
            gp(abs(gp)<10*eps) = 0; 
            gp(ceil(end/2)+1:end) = -flipud(gp(1:floor(end/2)));
            w(ceil(end/2)+1:end) = flipud(w(1:floor(end/2)));

        end  
    end
end
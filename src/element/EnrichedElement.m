%% EnrichedElement class
%
% This class defines a enriched finite element
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
% Initially prepared for the course CIV 2801 - Fundamentos de Computação
% Gráfica, 2022, second term, Department of Civil Engineering, PUC-Rio.
%
%% Class definition
classdef EnrichedElement < RegularElement
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        idEnr              = [];            % Vector identifying the intersections
        fracture           = [];            % Object fracture
        glw                = [];            % Vector of the enhancement degrees of freedom
        nglw               = 0;             % Number of enhancement dof
        beta               = 0.5;           % Ratio to compute the reference point (default is the mid point)
        enrVar             = 'w';           % Enrichment variable: 'w' or 'alpha'
        subDivInt          = false;         % Flag to apply a subdivision of the domain to define the quadrature points
        stretch            = false;         % Flag to indicate if the stretch part of the mapping matrix will be considered
        jumpOrder          = 1;             % Order of the interpolation of the jump displacement field
        staticCondensation = false          % Flag to apply a static condensation of the additional dofs
    end
    %% Constructor method
    methods
        function this = EnrichedElement(type, node, elem, anm, t, matModel, mat, intOrder, gla, fracture, glw, subDivInt, stretch, enrVar, jumpOrder, staticCondensation)
            this = this@RegularElement(type, node, elem, anm, t, matModel, mat, intOrder, gla);
            if nargin > 6
                this.fracture           = fracture;
                this.glw                = glw;
                if staticCondensation == true
                    this.gle = gla;
                else
                    this.gle = [gla glw];
                end
                this.ngla               = length(this.gla);
                this.nglw               = length(this.glw);
                this.ngle               = length(this.gle);
                this.subDivInt          = subDivInt;
                this.stretch            = stretch;
                this.enrVar             = enrVar;
                this.jumpOrder          = jumpOrder;
                this.staticCondensation = staticCondensation;
                this.createResults();
            end
        end
    end
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)

        % -----------------------------------------------------------------
        % Compute the matrix Gr. This matrix discretizes the bounded part
        % of the enhanced strains wrt to the vector of enhanced dof (w).
        [Gr, Gv] = enhancedStrainCompatibilityMtrcs(this, B, Xn)

    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % This function assembles the element stiffness matrix and internal
        % force vector
        % 
        % Input:
        %   dUe: vector with increment of the nodal displacement vector
        %        associated with the element. The first ngla components of
        %        this vector are associated with the regular dofs and the
        %        final components are associated to the enrichment dofs.
        %
        % Output:
        %   ke : element stiffness matrix
        %   fe : element internal force vector
        %
        function [ke,fe] = elementKeFint(this, dUe)

            % Initialize the element's stiffness matrix
            kaa = zeros(this.ngla, this.ngla);
            kaw = zeros(this.ngla, this.nglw);
            kwa = zeros(this.nglw, this.ngla);
            kww = zeros(this.nglw, this.nglw);

            % Initialize the internal force sub-vectors
            fa = zeros(this.ngla, 1);
            fw = zeros(this.nglw, 1);

            % Compute the increment of the enrichment dofs
            dwe = this.computeIncrEnrichmentDofs(dUe);

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints
            
                % Compute the B matrix at the int. point and the detJ
                [B, detJ] = this.shape.BMatrix(this.node,this.intPoint(i).X);

                % Compute the matrices Gr and Gv
                [Gr, Gv] = this.enhancedStrainCompatibilityMtrcs(B,this.intPoint(i).X);
        
                % Compute the increment of the strain vector
                dStrain = B*dUe(1:this.ngla) + Gr*dwe;
        
                % Compute the stress vector and the constitutive matrix
                [stress,D] = this.intPoint(i).constitutiveModel(dStrain);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
        
                % Numerical integration of the stiffness sub-matrices
                kaa = kaa + B' * D * B  * c;
                kaw = kaw + B' * D * Gr * c;
                kwa = kwa + Gv'* D * B  * c;
                kww = kww + Gv'* D * Gr * c;

                % Numerical integration of the internal force sub-vectors
                fa = fa + B'  * stress * c;
                fw = fw + Gv' * stress * c;
            end

            % Get the discontinuity stiffness matrix and internal force
            % vector
            [kd, fd] = this.fracture.elementKeFint(dwe,this.enrVar);

            % Add the fracture stiffness contribution
            kww = kww + kd;
            fw  = fw  + fd;

            % Assemble the element stiffness matrix and internal force
            % vector
            [ke,fe] = this.assembleElemKeFe(kaa,kaw,kwa,kww,fa,fw);
            
        end

        % -----------------------------------------------------------------
        % Function to assemble the elements stiffness matrix and internal
        % force vector
        function [ke,fe] = assembleElemKeFe(this,kaa,kaw,kwa,kww,fa,fw)

            if this.staticCondensation == true

                ke = kaa - kaw*(kww\kwa);
                fe = fa  - kaw*(kww\fw);

            else

                ke = [kaa, kaw;
                      kwa, kww]; 
                fe = [fa; fw];

            end
        end

        % -----------------------------------------------------------------
        % Function to compute the increment of the enrichment dofs
        function dWe = computeIncrEnrichmentDofs(this,dUe)

            if this.staticCondensation == true

                dWe = this.staticCondensation_ComputeIncrEnrDofs(dUe);

            else

                dWe = dUe((1+this.ngla):end);

            end

        end

        %------------------------------------------------------------------
        % Compute the increment of the enrichment dof based on the static
        % condensation process.
        function dwe = staticCondensation_ComputeIncrEnrDofs(this,dae)

            % Initialize the element's stiffness matrix
            kwa = zeros(this.nglw, this.ngla);
            kww = zeros(this.nglw, this.nglw);

            % Initialize the internal force sub-vectors
            fw = zeros(this.nglw, 1);

            % To compute the necessary stiffness matrices and internal
            % force vector, initialize the strain vector as zero.
            dwe_0     = zeros(this.nglw,1);
            dStrain_0 = zeros(3,1);

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints
            
                % Compute the B matrix at the int. point and the detJ
                [B, detJ] = this.shape.BMatrix(this.node,this.intPoint(i).X);

                % Compute the matrices Gr and Gv
                [Gr, Gv] = this.enhancedStrainCompatibilityMtrcs(B,this.intPoint(i).X);
        
                % Compute the stress vector and the constitutive matrix
                [stress,D] = this.intPoint(i).constitutiveModel(dStrain_0);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
        
                % Numerical integration of the stiffness sub-matrices
                kwa = kwa + Gv'* D * B  * c;
                kww = kww + Gv'* D * Gr * c;

                % Numerical integration of the internal force sub-vectors
                fw = fw + Gv' * stress * c;
            end

            % Get the discontinuity stiffness matrix and internal force
            % vector
            [kd, fd] = this.fracture.elementKeFint(dwe_0,this.enrVar);

            % Add the fracture stiffness contribution
            kww = kww + kd;
            fw  = fw  + fd;

            % Assemble the element stiffness matrix and internal force
            % vector
            dwe = kww\(fw - kwa*dae);

        end

        % -----------------------------------------------------------------
        % Compute the matrix that discretizes the bounded enhanced strain
        % field with the enrichment dofs.
        function Gkin = enhancedKinematicCompatibilityMtrx(this, B, Xn)

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
            Gkin = B*(h*I - Hd)*M;

        end

        % -----------------------------------------------------------------
        % Compute the matrix Gv. This matrix discretized the bounded part
        % of the virtual enhanced strains wrt to the vector of enhanced
        % dof (w). 
        function Gv = enhancedStaticCompatibilityMtrx(this, ~, Xn)

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
            Se = this.fracture.transformAlphaToW();
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
            
            % Numerical integration of the stiffness matrix components
            Gramm = zeros(3,3);
            for i = 1:this.nIntPoints

                % Integration point in the global cartesian coordinate
                % system
                [~, detJ] = this.shape.BMatrix(this.node,this.intPoint(i).X);
                X = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Integration point in the element local cartesian
                % coordinate system, with origin located at its
                % centroid
                Xrel = X - X0;
        
                % Compute the elastic constitutive matrix
                dA = [  1.0          Xrel(1)           Xrel(2);
                       Xrel(1)    Xrel(1)*Xrel(1)   Xrel(2)*Xrel(1);
                       Xrel(2)    Xrel(1)*Xrel(2)   Xrel(2)*Xrel(2)];
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ;
        
                % Numerical integration of the stiffness matrix Kaa
                Gramm = Gramm + dA * c;
            
            end

            % Numerical integration
            vec = zeros(3,1);
            for i = 1:this.fracture.nIntPoints
                WdetJ = this.fracture.intPoint(i).w * ld;
                Nw    = this.fracture.shape.shapeFncMtrx(this.fracture.intPoint(i).X);
                X     = Nw*[this.fracture.node(1,:),this.fracture.node(2,:)]';
                Xrel  = X - X0;
                s     = this.fracture.m*(X - this.fracture.Xref');
                dVec  = [s^k; (s^k)*Xrel(1); (s^k)*Xrel(2)];
                vec   = vec + dVec * WdetJ;
            end

            % Compute the coefficients
            C = Gramm \ vec;

        end



        % -----------------------------------------------------------------
        % Mapping matrix associated to a element. This matrix is
        % constructed by stacking by rows the mapping matrices evaluated at
        % the element's nodes.
        function Me = elementMappingMtrx(this)
            
            % Initialize the element's mapping matrix
            Me = zeros(this.nnd_el*this.ndof_nd,this.nglw);

            % Get the Poisson ratio of the material
            nu = this.mat(2);

            % Construct the matrix by stacking by-rows the mapping matrix
            % evaluated at each node
            for i = 1:this.nnd_el

                % Compute the mapping matrix at node X
                X  = this.node(i,:);
                Mi = this.fracture.jumpTransmissionMtrx(X,this.enrVar,this.stretch,nu);

                % Assemble the element mapping matrix
                Me((2*i-1):(2*i),:) = Mi;
            end

        end

        % -----------------------------------------------------------------
        % Heaviside function evaluated at the point (or set of points) X,
        % wrt to the reference point of the fracture (Xref).
        function h = heavisideFnc(this,X)

            % Reference point at the fracture that crosses the element
            Xref = this.fracture.Xref;

            % Fracture normal vector
            n = this.fracture.n;

            % Relative distance of the nodes of the element wrt to the reference point
            DX = X - repmat(Xref,size(X,1),1);
            
            % Heaviside function evaluated in the nodes of the element
            h = max(sign(DX*n'),0.0);

        end

        % -----------------------------------------------------------------
        % Diagonal matrix with the Heaviside function evaluated at the
        % nodes associeted to each regular dof of the element, wrt to the
        % reference point of the fracture (Xref).
        function Hde = heavisideMtrx(this)
 
            % Heaviside function evaluated in the nodes of the element
            h = this.heavisideFnc(this.node);

            % Create the vector
            hv = repmat(h,[1,this.ndof_nd])';
            
            % Matrix with the Heaviside function evaluated in the node of the dofs
            Hde = diag(reshape(hv,numel(hv),1));

        end

        %------------------------------------------------------------------
        % Function to compute the displacement field inside a given ISOQ4 
        % element.
        function u = displacementField(this,X)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   u   : displacement vector evaluated in "X"
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Heaviside function evaluated in X
            h = this.heavisideFnc(X);

            % Compute the Heaviside matrix
            Hd = this.heavisideMtrx();

            % Compute the mapping matrix
            M = this.elementMappingMtrx();

            % Get enrichment dofs
            
            % Displacement field
            u = Nm*this.ue(1:this.ngla) + Nm*(h*eye(size(Hd)) - Hd)*M*this.ue(this.ngla+1:end);
        
        end

        %------------------------------------------------------------------
        % Function to get the nodes that will be used to plot the results.
        function [resNodes,fractNodes] = getResultNodes(this)

            % Define the edges that the fracture crosses. The following
            % methodology works for one fracture only, then to take into
            % account multiple fractures, a loop needs to be performed.
            
            % Order the nodes in counterclockwise order
            Nodes = [this.node; this.fracture.node];
            order = this.sortCounterClockWise(Nodes);

            % Get the position of the fracture nodes in the order
            nf1 = find(order == size(this.node,1)+1);
            nf2 = find(order == size(this.node,1)+2);

            % The nodes before and after each node of the fracture define
            % the edge of the element.
            % Edge associated to the first node:
            if nf1 == 1
                edge1 = [order(end) order(nf1+1)];
            elseif nf1 == length(order)
                edge1 = [order(nf1-1) order(1)];
            else
                edge1 = [order(nf1-1) order(nf1+1)];
            end
            % Edge associated to the second node:
            if nf2 == 1
                edge2 = [order(end) order(nf2+1)];
            elseif nf2 == length(order)
                edge2 = [order(nf2-1) order(1)];
            else
                edge2 = [order(nf2-1) order(nf2+1)];
            end

            % Length of each edge:
            len1 = sqrt((Nodes(edge1(2),1)-Nodes(edge1(1),1))^2+(Nodes(edge1(2),2)-Nodes(edge1(1),2))^2);
            len2 = sqrt((Nodes(edge2(2),1)-Nodes(edge2(1),1))^2+(Nodes(edge2(2),2)-Nodes(edge2(1),2))^2);

            % Compute the vectors with the origin at the nodes of the
            % fracture oriented along the edges.
            vf1 = Nodes(edge1(1),:) - this.fracture.node(1,:);
            vf2 = Nodes(edge2(1),:) - this.fracture.node(2,:);
            vf1 = vf1/norm(vf1);
            vf2 = vf2/norm(vf2);

            % Compute the increment vector:
            df1 = (len1/1000)*dot(this.fracture.n,vf1)*vf1;
            df2 = (len2/1000)*dot(this.fracture.n,vf2)*vf2;

            % Find a node along the Omega^plus region along each edge:
            n1_plus = this.fracture.node(1,:) + df1;
            n2_plus = this.fracture.node(2,:) + df2;

            % Find a node along the Omega^minus region along each edge:
            n1_minus = this.fracture.node(1,:) - df1;
            n2_minus = this.fracture.node(2,:) - df2; 
           
            % Get the order of the nodes for the element connectivity be
            % done in a counterclockwise way
            fractNodes = [n1_plus;n1_minus; n2_plus;n2_minus];
            resNodes = [this.node; fractNodes];
            orderResNodes = this.sortCounterClockWise(resNodes);
            resNodes = resNodes(orderResNodes,:);

            % Nodes for the fracture result object
            orderFract = this.sortCounterClockWise(fractNodes);
            fractNodes = fractNodes(orderFract,:);

        end

        %------------------------------------------------------------------
        % Function to create an array of result objects based on the
        % children elements.
        function createResults(this)
            [resNodes,fractNodes] = this.getResultNodes();
            nNodes      = size(resNodes,1);
            nFractNodes = size(fractNodes,1);
            this.result = Result(resNodes ,1:nNodes ,0.0*ones(nNodes,1) ,'');
            this.fracture.result = Result(fractNodes ,1:nFractNodes ,0.0*ones(nNodes,1) ,'');
        end

    end

    methods(Static)

        %------------------------------------------------------------------
        % This function sorts counterclockwise a set of nodes.
        % It uses as a reference point the centroid defined by the nodes.
        function order = sortCounterClockWise(NODE)
            
            % Centroid coordinate
            cx = mean(NODE(:,1));
            cy = mean(NODE(:,2));
            
            % Compute the angle that the relative vector of the vertices 
            % from the centroid has with the horizontal axis
            a = atan2(NODE(:,2) - cy, NODE(:,1) - cx);
            
            % Sort the angles
            [~, order] = sort(a);
            
        end

    end
end
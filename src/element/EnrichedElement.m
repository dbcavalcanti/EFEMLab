%% EnrichedElement class
%
% This class defines a enriched finite element (consider a ISOQ4 element)
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
        idEnr     = [];            % Vector identifying the intersections
        fracture  = [];            % Object fracture
        glw       = [];            % Vector of the enhancement degrees of freedom
        ngla      = 0;             % Number of regular dof
        nglw      = 0;             % Number of enhancement dof
        ngle      = 0;             % Number of total dof
        beta      = 0.5;           % Ratio to compute the reference point (default is the mid point)
        subDivInt = false;         % Flag to apply a subdivision of the domain to define the quadrature points
        stretch   = false;         % Flag to indicate if the stretch part of the mapping matrix will be considered
        jumpOrder = 1;             % Order of the interpolation of the jump displacement field
    end
    %% Constructor method
    methods
        function this = EnrichedElement(type, node, elem, anm, t, matModel, mat, intOrder, gla, fracture, glw, subDivInt, stretch, jumpOrder)
            this = this@RegularElement(type, node, elem, anm, t, matModel, mat, intOrder, gla);
            if nargin > 6
                this.fracture  = fracture;
                this.glw       = glw;
                this.gle       = [gla glw];
                this.ngla      = length(this.gla);
                this.nglw      = length(this.glw);
                this.ngle      = length(this.gle);
                this.subDivInt = subDivInt;
                this.stretch   = stretch;
                this.jumpOrder = jumpOrder;
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
        Gr = enhancedStrainCompatibilityMtrx(this, B, Xn)

        % -----------------------------------------------------------------
        % Compute the matrix Gv. This matrix discretizes the bounded part
        % of the virtual enhanced strains wrt to the vector of enhanced
        % dof (w). 
        Gv = enhancedStressCompatibilityMtrx(this, B, Xn)

    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % This function assembles the element's stiffness matrix
        function ke = elementStiffnessMtrx(this, dUe)

            % Compute the stiffness sub-matrices
            [kaa, kaw, kwa, kww] = this.computeElemStiffnessSubMatrices(dUe);

            % Assemble the element stiffness matrix
            ke = [kaa, kaw;
                  kwa, (kww+this.fracture.kd)];
            
        end

        %------------------------------------------------------------------
        % This function computes the element's stiffness sub-matrices
        function [kaa, kaw, kwa, kww] = computeElemStiffnessSubMatrices(this,dUe)

            % Initialize the element's stiffness matrix
            kaa = zeros(this.ngla, this.ngla);
            kaw = zeros(this.ngla, this.nglw);
            kwa = zeros(this.nglw, this.ngla);
            kww = zeros(this.nglw, this.nglw);
             
            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints
            
                % Compute the B matrix at the int. point and the detJ
                [B, detJ] = this.shape.BMatrix(this.node,this.intPoint(i).X);

                % Compute the matrix Gr
                Gr = this.enhancedStrainCompatibilityMtrx(B,this.intPoint(i).X);

                % Compute the matrix Gv
                Gv = this.enhancedStressCompatibilityMtrx(B,this.intPoint(i).X);
        
                % Compute the increment of the strain vector
                dStrain = B*dUe(1:this.ngla);
        
                % Compute the elastic constitutive matrix
                De = this.intPoint(i).getConstitutiveMtrx(dStrain);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
        
                % Numerical integration of the stiffness sub-matrices
                kaa = kaa + B' * De * B  * c;
                kaw = kaw + B' * De * Gr * c;
                kwa = kwa + Gv'* De * B  * c;
                kww = kww + Gv'* De * Gr * c;
            end

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the translation part of the jump displacement
        % field
        function Mtr = mappingTranslation(this)

            % Fracture geometric properties
            m    = this.fracture.m;
            ld   = this.fracture.ld;
            Xref = this.fracture.Xref;
            X1   = this.fracture.node(1,:);
            X2   = this.fracture.node(2,:);

            % Tangential relative coordinate
            s1 = m*(X1' - Xref');
            s2 = m*(X2' - Xref');

            Mtr = [ s2     0.0   -s1    0.0;
                    0.0     s2    0.0  -s1]/ld; 

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the relative rotation part of the jump 
        % displacement field. Evaluated at one point X.
        function Mrot = mappingRelativeRotation(this,X)
            
            % Fracture geometric properties
            m    = this.fracture.m;
            ld   = this.fracture.ld;
            Xref = this.fracture.Xref;

            % Mapping matrix
            Mrot = [-(X(2)-Xref(2))*m(2)   (X(1)-Xref(1))*m(2);
                     (X(2)-Xref(2))*m(1)  -(X(1)-Xref(1))*m(1);
                     (X(2)-Xref(2))*m(2)  -(X(1)-Xref(1))*m(2);
                    -(X(2)-Xref(2))*m(1)   (X(1)-Xref(1))*m(1)]' / ld;

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the stretching part of the jump 
        % displacement field. Evaluated at one point X.
        function Ms = mappingStretching(this,X)
            
            % Fracture geometric properties
            m    = this.fracture.m;
            ld   = this.fracture.ld;
            Xref = this.fracture.Xref;
            cs   = m(1);
            sn   = m(2);

            % Tangential relative coordinate
            s = m*(X' - Xref');

            % Mapping matrix
            Ms = [-cs*cs   -sn*cs;
                  -sn*cs   -sn*sn;
                   cs*cs    sn*cs;
                   sn*cs    sn*sn]' *(s/ld);

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the rigid-body movement due to the rotation.
        % Evaluated at one point X
        function Mrb = mappingRigidBody(this,X)
            
            % Mapping due to the translation
            Mrb  = this.mappingTranslation();

            % Add the contribution of the relative rotation
            if this.jumpOrder == 1
                Mrb = Mrb + this.mappingRelativeRotation(X);
            end

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the non-rigid body movement due to the rotation.
        % Evaluated at one point X
        function Mnrb = mappingNonRigidBody(this,X)

            % Mapping due to the stretching
            Mnrb = this.mappingStretching(X);

        end

        % -----------------------------------------------------------------
        % Mapping matrix associated to a element. This matrix is
        % constructed by stacking by rows the mapping matrices evaluated at
        % the element's nodes.
        function Me = elementMappingMtrx(this)
            
            % Initialize the element's mapping matrix
            Me = zeros(this.nnd_el*this.ndof_nd,4);

            % Construct the matrix by stacking by-rows the mapping matrix
            % evaluated at each node
            for i = 1:this.nnd_el

                % Compute the mapping matrix at node X
                X  = this.node(i,:);
                Mi = this.mappingRigidBody(X);
                if (this.jumpOrder > 0) && (this.stretch == true)
                    Mi = Mi + this.mappingNonRigidBody(X);
                end

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
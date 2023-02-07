%% Element class
%
% This class defines a finite element model (consider a ISOQ4 element)
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
classdef RegularElement < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        type       = 'ISOQ4';      % type of element
        shape      = [];            % Object of the Shape class
        node       = [];            % Nodes of the fem mesh
        connect    = [];            % Nodes connectivity
        anm        = 'PlaneStress'; % Analysis model
        t          = 1.0;           % Thickness
        matModel   = 'elastic';     % Material model
        mat        = [];            % Vector with material properties
        intOrder   = 2;             % Order of the numerical integration
        nnd_el     = 4;             % Number of nodes per element
        ndof_nd    = 2;             % Number of dof per node
        gla        = [];            % Vector of the regular degrees of freedom
        gle        = [];            % Vector of the degrees of freedom
        ngla       = 0;             % Number of regular dof
        ngle       = 0;             % Number of total dof
        ue         = [];            % Element's displacement vector
        nIntPoints = 1;             % Number of integration points
        intPoint   = [];            % Vector with integration point objects
        result     = [];            % Result object to plot the results
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement(type, node, elem, anm, t, ...
                matModel, mat, intOrder, gla)
            if (nargin > 0)
                if strcmp(type,'ISOQ4')
                    this.shape = Shape_ISOQ4();
                elseif strcmp(type,'CST')
                    this.shape = Shape_CST();
                end
                this.node     = node;
                this.nnd_el   = size(node,1);
                this.connect  = elem;
                this.t        = t;
                this.matModel = matModel;
                this.mat      = mat;
                this.anm      = anm;
                this.intOrder = intOrder;
                this.gla      = gla;
                this.ngla     = length(this.gla);
                this.gle      = gla;
                this.ngle     = length(this.gle);
                this.result   = Result(this.node,1:length(this.connect),0.0*ones(this.nnd_el,1),'Model');
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(this.intOrder,this);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                if strcmp(this.matModel,'elastic')
                    constModel = Material_Elastic(this.mat, this.anm);
                end
                intPts(i) = IntPoint(X(:,i),w(i),this.anm, constModel);
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % This function assembles the element stiffness matrix and internal
        % force vector
        % 
        % Input:
        %   dUe: vector with increment of the nodal displacement vector
        %        associated with the element
        %
        % Output:
        %   ke : element stiffness matrix
        %   fe : element internal force vector
        %
        function [ke,fe] = elementKeFint(this,dUe)

            % Initialize the element stiffness matrix and internal force
            % vector
            ke = zeros(this.ndof_nd*this.nnd_el);
            fe = zeros(this.ndof_nd*this.nnd_el,1);
            
            % Numerical integration of the stiffness matrix and internal
            % force vector
            for i = 1:this.nIntPoints
            
                % Compute the B-matrix and detJ at the integration point 
                [B,detJ] = this.shape.BMatrix(this.node,this.intPoint(i).X);

                % Compute the increment of the strain vector
                dStrain = B*dUe;
        
                % Compute the stress vector and the constitutive matrix
                [stress,D] = this.intPoint(i).constitutiveModel(dStrain);
        
                % Numerical integration term
                c = this.intPoint(i).w * detJ * this.t;
        
                % Numerical integration of the stiffness matrix and the
                % internal force vector
                ke = ke + B' * D * B  * c;
                fe = fe + B' * stress * c;

            end
            
        end

        %------------------------------------------------------------------
        % Function to compute the displacement field in the element.
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
            
            % Regular displacement field
            u = Nm*this.ue;
        
        end

        %------------------------------------------------------------------
        % Function to compute the area of the element domain
        function A = getDomainArea(this)
            nd = [this.connect(2) this.connect(1)];
            A =-(this.node(nd(2),2)+this.node(nd(1),2))*(this.node(nd(2),1)-this.node(nd(1),1))/2;
            for j = 2:length(this.connect)
                nd = [this.connect(j-1) this.connect(j)];
                A  = A - (this.node(nd(2),2)+this.node(nd(1),2))*...
                    (this.node(nd(2),1)-this.node(nd(1),1))/2;
            end

        end
        

    end

    %% Public methods associated with the pos-processing
    methods

        %------------------------------------------------------------------
        % Update the result's object vertices property
        % If the 'Undeformed' configuration is selected, nothing needs to
        % be done.
        function updateResultVertices(this,configuration)
            if strcmp(configuration,'Deformed')
                Nodes = this.getDeformedConfiguration();
                this.result.setVertices(Nodes);
            end  
        end

        %------------------------------------------------------------------
        % Update the result's object vertices property
        % If the 'Undeformed' configuration is selected, nothing needs to
        % be done.
        function updateResultFaces(this,faces)
            this.result.setFaces(faces);
        end
        
        %------------------------------------------------------------------
        % Update the result's object vertices data property
        function updateResultVertexData(this,type)
            this.result.setDataLabel(type);
            switch type
                case 'Ux'
                    ndResults = this.getNodalDisplacementField(1);
                case 'Uy'
                    ndResults = this.getNodalDisplacementField(2);
                case 'Sxx'
                    ndResults = this.getNodalStressField(1);
                case 'Syy'
                    ndResults = this.getNodalStressField(2);
                case 'Sxy'
                    ndResults = this.getNodalStressField(3);
            end
            this.result.setVertexData(ndResults);
        end

        %------------------------------------------------------------------
        % Get the specified displacement field component at the nodes of
        % the model.
        function NodeDef = getDeformedConfiguration(this)
            NodeDef = this.node;
            Ue      = reshape(this.ue,[size(this.node,2),size(this.node,1)])';
            for i = 1:this.nnd_el
                NodeDef(i,:) = NodeDef(i,:) + Ue(i,:);
            end
        end

        %------------------------------------------------------------------
        % Get the specified displacement field component at the nodes of
        % the model.
        function Ui = getNodalDisplacementField(this,component)
            Ui = this.ue(component:this.ndof_nd:this.nnd_el*this.ndof_nd);
        end

    end

end
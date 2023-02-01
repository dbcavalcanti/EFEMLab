%% Fracture class
%
% This class defines a fracture
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
classdef Fracture < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        shape      = [];            % Object of the Shape class
        idelem     = 0;             % Identify the element that the fracture belongs
        node       = [];            % Nodes of the fem mesh
        connect    = [];            % Nodes connectivity
        m          = [];            % Tangent orientation vector
        n          = [];            % Normal orientation vector
        ld         = 0.0;           % Fracture length
        Xref       = [];            % Reference point
        t          = 1.0;           % Thickness
        matModel   = 'elastic';     % Material model
        mat        = [];            % Vector with material properties
        penal      = false;         % Flag for applying a penalization on compression
        nnd_el     = 2;             % Number of nodes per element
        ndof_nd    = 2;             % Number of dof per node
        ndof       = 1;             % Number of dofs
        glw        = [];            % Vector of the degrees of freedom
        nIntPoints = 2;             % Number of integration points
        intPoint   = [];            % Vector with integration point objects
        result     = [];            % Result object
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Fracture(node, elem, t, matModel, mat, glw, penal)
            if (nargin > 0)

                this.node     = node;
                this.connect  = elem;
                this.t        = t;
                this.matModel = matModel;
                this.mat      = mat;
                this.glw      = glw;
                this.penal    = penal;
                this.ndof     = length(glw);
                this.shape    = Shape_Bar();

                % Initialize the geometry properties
                this.initializeGeometry();

                % Initialize the integration points
                this.initializeIntPoints();
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the fracture element geometry. Computes the length,
        % the orientation vectors and the reference point.
        function initializeGeometry(this)

            % Fracture length
            dx = this.node(2,1) - this.node(1,1);
            dy = this.node(2,2) - this.node(1,2);
            this.ld = sqrt(dx.^2 + dy.^2);
            
            % Fracture orientation
            sn = dy./this.ld;
            cs = dx./this.ld;
            
            % Tangential vector to the discontinuity
            this.m = [ cs   sn];  
            
            % Normal vector to the discontinuity
            % Defined considering n = ez x m, where ez = [0 0 1]
            this.n = [-sn   cs];

            % Reference point
            this.Xref = 0.5*(this.node(1,:) + this.node(2,:));
            
        end

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints();

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints

                % Initialize the constitutive model
                if strcmp(this.matModel,'elastic')
                    constModel = MaterialInterface_Elastic(this.mat, this.penal);
                elseif strcmp(this.matModel,'isotropicDamage')
                    constModel = MaterialInterface_IsotropicDamage(this.mat, this.penal);
                end

                % Create the integration points
                intPts(i) = IntPoint(X(:,i),w(i),'Interface', constModel);

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
            ke = zeros(this.ndof,this.ndof);
            fe = zeros(this.ndof,1);

            % Compute the rotation matrix
            R = this.rotationMtrx();

            % Rotate the jump displacement vector to the local system
            dUe = R * dUe;

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints

                % Shape function matrix
                Nw = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Evaluate the jump at the integration point in the local
                % coordinate system
                dw = Nw * dUe;
           
                % Compute the stress vector and the constitutive matrix
                [td,T] = this.intPoint(i).constitutiveModel(dw);
        
                % Numerical integration term
                c = this.intPoint(i).w * this.ld * this.t;
        
                % Numerical integration of the stiffness matrix and the
                % internal force vector
                ke = ke + Nw' * T  * Nw * c;
                fe = fe + Nw' * td * c;

            end

            % Rotate to the global coordinate system
            ke = R' * ke * R;
            fe = R' * fe;
            
        end

        %------------------------------------------------------------------
        % This function computes the element's rotation matrix. Change from
        % the local coordinate system mn to the global system xy
        function R = rotationMtrx(this)

            % Rotation of a point
            r = [ this.m(1,1)   this.m(1,2);
                 -this.m(1,2)   this.m(1,1) ];

            % Rotation matrix of the element (2 points)
            R = blkdiag(r,r);

        end

    end

end
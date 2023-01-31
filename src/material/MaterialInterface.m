%% TractionConstitutiveLaw class
%
% This class defines an abstract traction-displacement constitutive law
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
classdef MaterialInterface < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        kn = 0.0;   
        ks = 0.0;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface(node, elem, t, mat, glw)
            if (nargin > 0)
                this.node    = node;
                this.connect = elem;
                this.t       = t;
                this.mat     = mat;
                this.glw     = glw;
                % Initialize the geometry properties of the fracture
                % element
                this.initializeGeometry();
                % Initialize the fracture stiffness matrix
                this.stiffnessMtrx();
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the fracture element geometry
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
            this.n = [-sn   cs];

            % Reference point
            this.Xref = 0.5*(this.node(1,:) + this.node(2,:));
            
        end

        %------------------------------------------------------------------
        % This function computes the fracture's stiffness matrix in the 
        % local system
        function localStiffnessMtrx(this)

            % Initialize the element's stiffness matrix
            this.kd = zeros(this.ndof_nd*this.nnd_el);
            
            % Compute the integration points
            [w,x] = this.getlineQuadrature_GaussLobatto();
            
            % Numerical integration of the stiffness matrix components
            Td = this.constitutiveMtrx();
            for i = 1:length(x)
                WdetJ    = w(i)*this.ld;
                Nw       = this.shapeFncFracture(x(i));
                this.kd  = this.kd + Nw' * Td * Nw * WdetJ * this.t;
            end

        end

        %------------------------------------------------------------------
        % This function computes the element's stiffness matrix 
        function stiffnessMtrx(this)

            % Compute the element's stiffness matrix in the local system
            this.localStiffnessMtrx();

            % Compute the rotation matrix
            R = this.rotationMtrx();
            
            % Rotate the stiffness matrix to the global coordinate system
            this.kd = R' * this.kd * R;

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

        %------------------------------------------------------------------
        % This function return the global cartesian coordinates of a point inside a -
        % isoparametric linear triangular element for a given point in the natural 
        % coordinate system.
        function Td = constitutiveMtrx(this)
        
            % Elastic material properties
            ks = this.mat(1);
            kn = this.mat(2);
            
            % CCnstitutive matrix
            Td = [ ks   0.0;
                  0.0   kn ];
        end

    end

    %% Public static methods
    methods (Static)

        %------------------------------------------------------------------
        % Get Gauss-Lobatto quadrature of order 2
        function [w,gp] = getlineQuadrature_GaussLobatto()

            gp = [0.0   1.0];   % Weights
            w  = [0.5   0.5];   % Coordinates in the natural system    

        end
    
        %------------------------------------------------------------------
        % Returns the shape function of a linear one-dimensional
        % isoparametricelement.
        function N = shapeFncFracture(s)

            N = [(1.0-s) ,    0.0  ,  s  , 0.0;
                    0.0  , (1.0-s) , 0.0 ,  s];
            
        end

    end
end
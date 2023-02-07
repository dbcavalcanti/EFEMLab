%% Model class
%
% This class defines a finite element model that has a strong discontinuity
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%
% Updates:
%   January 2023:
%      Possibility to choose between the Discrete Strong Discontinuiy 
%   Approach (DSDA) and the Generalized Strong Discontinuity Approach 
%   (GSDA) by choosing to consider or not the non-rigid body part of the
%   mapping matrix.
%      Possibility to choose between the Kinematically Optimal Symmetric
%   Formulation (KOS) and the Kinematically and Statically Optimal
%   Non-Symmetric Formulation (KSON).
% 
%
%% Class definition
classdef Model < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        NODE                = [];            % Nodes of the fem mesh
        ELEM                = [];            % Nodes connectivity
        t                   = 1.0;           % Thickness
        matModel            = 'elastic';     % Continuum constitutive law
        mat                 = [];            % Vector with material properties
        anm                 = 'PlaneStress'; % Analysis model identification
        type                = 'ISOQ4';       % Typf of element used
        SUPP                = [];            % Matrix with support conditions
        LOAD                = [];            % Matrix with load conditions
        PRESCDISPL          = [];            % With with prescribed displacements
        intOrder            = 2;             % Number of integration points
        NODE_D              = [];            % Nodes of the fractures
        FRACT               = [];            % Fractures' nodes connectivity
        tractionLaw         = 'elastic';     % Discontinuity traction separation law
        tractionLawPenal    = false;         % Flag for a penalination of compression in the traction separation law
        matfract            = [];            % Vector with the fracture cohesive parameters
        IDenr               = [];            % Matrix identifying the intersections
        enhancementType     = 'KOS';         % String with the type of the SDA formulation ('KOS', 'KSON')
        subDivInt           = false;         % Flag for applying a sub-division of the domain to perform the numerical integration
        stretch             = false;         % Flag for considering the stretch part of the mapping matrix
        enrVar              = 'w';           % Chosen enrichment variable: 'w' or 'alpha'
        lvlEnrVar           = 'global';   % Level of the enrichment dofs ('local or 'global')
        jumpOrder           = 1;             % Order of the interpolation of the jump displacement field
        staticCondensation  = false;         % Flag for applying a static condensation of the additional dofs
        nnodes              = 1;             % Number of nodes
        nfracnodes          = 0;             % Number of fracture nodes
        nelem               = 1;             % Number of elements
        nnd_el              = 4;             % Number of nodes per element
        ndof_nd             = 2;             % Number of dof per node
        ndof_frac           = 4;             % Number of dof per frac
        ndof                = 1;             % Number of regular degrees of freedom
        nTotDofs            = 0;             % Total number of degrees of freedom
        ndoffree            = 0;             % Number of free degrees of freedom
        ndoffixed           = 0;             % Number of fixed degrees of freedom
        enrDof              = [];            % Vector with all enrichment dofs
        enrFreeDof          = [];            % Vector with the free enrichment dofs
        totFreeDof          = [];            % Vector with the total free dofs
        ID                  = [];            % Each line of the ID matrix contains the global numbers for the node DOFs (DX, DY)
        IDfrac              = [];            % Each line of the ID matrix contains the global numbers for the node of the fracture DOFs (DX, DY)
        GLA                 = [];            % Matrix with the regular dof of each element
        GLW                 = [];            % Cell with the enhacement dof of each element
        F                   = [];            % Global force vector
        U                   = [];            % Global displacement vector
        element             = [];            % Array with the element's objects
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Model(NODE, ELEM, NODE_D, FRACT, t, matModel, ...
                mat, tractionLaw, tractionLawPenal, matfract, anm, type, SUPP, LOAD, ...
                PRESCDISPL, intOrder,enhancementType, subDivInt, ...
                stretch, enrVar, jumpOrder, lvlEnrVar, staticCondensation, IDenr)
            if (nargin > 0)
                this.NODE               = NODE;
                this.ELEM               = ELEM;
                this.type               = type;
                this.t                  = t;
                this.matModel           = matModel;
                this.mat                = mat;
                this.anm                = anm;
                this.SUPP               = SUPP;
                this.LOAD               = LOAD;
                this.PRESCDISPL         = PRESCDISPL;
                this.intOrder           = intOrder;
                this.NODE_D             = NODE_D;
                this.FRACT              = FRACT;
                this.tractionLaw        = tractionLaw;
                this.tractionLawPenal   = tractionLawPenal;
                this.matfract           = matfract;
                this.enhancementType    = enhancementType;
                this.subDivInt          = subDivInt;
                this.stretch            = stretch;
                this.enrVar             = enrVar;
                this.lvlEnrVar          = lvlEnrVar;
                this.staticCondensation = staticCondensation;
                this.jumpOrder          = jumpOrder;
                this.IDenr              = IDenr;
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % this pre-computations
        function preComputations(this)
            
            % Initialize basic variables
            this.nnodes     = size(this.NODE,1);
            this.nfracnodes = size(this.NODE_D,1);
            this.nelem      = size(this.ELEM,1);     
            this.nnd_el     = size(this.ELEM,2);    
            this.ndof_nd    = 2;               
            this.ndof       = this.ndof_nd * this.nnodes; 

            % --- Assemble nodes DOF ids matrix ---------------------------
            %   Each line of the ID matrix contains the global numbers for 
            %   the node DOFs (DX, DY). Free DOFs are numbered first.
            
            % Initialize the ID matrix and the number of fixed dof
            this.ID = zeros(this.nnodes,this.ndof_nd);
            this.ndoffixed = 0;
            
            % Assemble the ID matrix
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    if (this.SUPP(i,j) == 1)
                        this.ndoffixed = this.ndoffixed + 1;
                        this.ID(i,j) = 1;
                    end
                end
            end
            
            % Number of free dof
            this.ndoffree = this.ndof - this.ndoffixed;
            
            % Initialize the counters
            countS = this.ndoffree;
            countF = 0;
            
            % Update the ID matrix with the free dof numbered first
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    if this.ID(i,j) == 0
                        countF  = countF + 1;
                        this.ID(i,j) = countF;
                    else
                        countS  = countS + 1;
                        this.ID(i,j) = countS;
                    end
                end
            end
            
            % Assemble the matrix with the degrees of freedom of each 
            % element
            this.GLA = zeros(this.nelem, this.nnd_el*this.ndof_nd);
            for el = 1:this.nelem
                this.GLA(el,:) = reshape(this.ID(this.ELEM(el,:),:)',1,...
                    this.nnd_el*this.ndof_nd);
            end

            % Initialize the cell with the element enrichment dofs
            % It is used a cell since there can be more than one fracture
            % embedded inside the element.
            % Initialize the matrix with the enhanced dof associated to
            % each node
            this.GLW = cell(1,this.nelem);
           
            if strcmp(this.lvlEnrVar,'global') && (this.jumpOrder == 1) && strcmp(this.enrVar,'w')
                
                this.IDfrac = zeros(size(this.NODE_D));
                count = this.ndof + 1;
                for i = 1:this.nfracnodes
                    this.IDfrac(i,:) = [count, count+1];
                    count = count + 2;
                end
                for el = 1:this.nelem
                    if sum(this.IDenr(el,:)) > 0
                        id = find(this.IDenr(el,:)==1);
                        for i = 1:length(id)
                            this.GLW{el} = [this.IDfrac(this.FRACT(id(i),1),:),...
                                            this.IDfrac(this.FRACT(id(i),2),:)];
                        end
                    end
                end
            
            elseif strcmp(this.lvlEnrVar,'local') && (this.jumpOrder == 1) && (strcmp(this.enrVar,'w') || strcmp(this.enrVar,'alpha'))

                count = this.ndof + 1;
                for el = 1:this.nelem
                    if sum(this.IDenr(el,:)) > 0
                        id = find(this.IDenr(el,:)==1);
                        for i = 1:length(id)
                            this.GLW{el} = [count count+1 count+2 count+3];
                            count = count + 4;
                        end
                    end
                end

            elseif strcmp(this.lvlEnrVar,'local') && (this.jumpOrder == 0) && (strcmp(this.enrVar,'w') || strcmp(this.enrVar,'alpha')) 

                count = this.ndof + 1;
                for el = 1:this.nelem
                    if sum(this.IDenr(el,:)) > 0
                        id = find(this.IDenr(el,:)==1);
                        for i = 1:length(id)
                            this.GLW{el} = [count count+1];
                            count = count + 2;
                        end
                    end
                end

            else

%                 fprintf(0,'Error. Check the selected EFEM formulation set up');
%                 error('EFEM formulation set up');

            end

            % Vector with all enrichment dofs
            this.enrDof     = unique(cell2mat(this.GLW))';
            this.enrFreeDof = this.enrDof;

            % Number of total dofs
            if this.staticCondensation == true
                this.nTotDofs = this.ndof;
                this.totFreeDof = 1:this.ndoffree;
            elseif this.staticCondensation == false
                this.nTotDofs = this.ndof + length(this.enrFreeDof);
                this.totFreeDof = [1:this.ndoffree, this.enrFreeDof'];
            end

            % Initialize the vector with the Element's objects
            elements(this.nelem,1) = Element(); 

            % Assemble the properties to the elements' objects
            for el = 1 : this.nelem
                if sum(this.IDenr(el,:)) == 0

                    elements(el) = RegularElement(...
                        this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                        this.anm, this.t, this.matModel, this.mat, this.intOrder,this.GLA(el,:));

                elseif sum(this.IDenr(el,:)) == 1

                    % Get which fractures are embedded in the element
                    id = find(this.IDenr(el,:)==1);

                    % Initialize the fracture object
                    if this.jumpOrder == 0
                        fract = Fracture_ConstantJump(this.NODE_D(this.FRACT(id,:),:),...
                            this.FRACT(id,:)+this.nnodes, this.t, this.tractionLaw, this.matfract(id,:), ...
                            this.GLW{el},this.tractionLawPenal);
                    elseif this.jumpOrder == 1
                        fract = Fracture_LinearJump(this.NODE_D(this.FRACT(id,:),:),...
                            this.FRACT(id,:)+this.nnodes, this.t, this.tractionLaw, this.matfract(id,:), ...
                            this.GLW{el},this.tractionLawPenal);
                    end

                    % Initialize the enriched elements using:
                    if strcmp(this.enhancementType,'KOS')

                        % The Kinematically Optimal Symmetric formulation
                        elements(el) = EnrichedElement_KOS(this.type,...
                            this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.anm,this.t,this.matModel, this.mat, this.intOrder,this.GLA(el,:),...
                            fract, this.GLW{el},this.subDivInt,...
                            this.stretch, this.enrVar,this.jumpOrder,this.staticCondensation);

                    elseif strcmp(this.enhancementType,'KSON')

                        % The Kinematically and Statically Optimal Nonsymmetric formulation
                        elements(el) = EnrichedElement_KSON(this.type,...
                            this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.anm,this.t, this.matModel,this.mat, this.intOrder,this.GLA(el,:),...
                            fract, this.GLW{el},this.subDivInt,...
                            this.stretch, this.enrVar,this.jumpOrder,this.staticCondensation);

                    end
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
            
            % Assemble load vector 
            this.F = zeros(this.ndof + length(this.enrDof),1);
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.F(this.ID(i,j)) = this.F(this.ID(i,j)) + ...
                        this.LOAD(i,j);
                end
            end
            
            % Initialize the displacement vector 
            this.U = zeros(this.ndof + length(this.enrDof),1);

            % Add the prescribed displacements
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.U(this.ID(i,j)) = this.U(this.ID(i,j)) + ...
                        this.PRESCDISPL(i,j);
                end
            end

        end

        %------------------------------------------------------------------
        % Global stiffness matrix and internal force vector
        function [K, Fint] = globalKF(this,dU)   
            
            % Initialize the global stiffness matrix
            K    = sparse(this.nTotDofs, this.nTotDofs);
            Fint = zeros(this.nTotDofs, 1);
            
            for el = 1:this.nelem

                % Get the vector of the element dof
                gle = this.element(el).type.gle;

                % Get the elements displacement vector
                dUe = dU(gle);
            
                % Get local stiffness matrix
                [ke,fe] = this.element(el).type.elementKeFint(dUe);
            
                % Assemble
                K(gle,gle) = K(gle,gle) + ke;
                Fint(gle) = Fint(gle) + fe;
                
            end
        end

      % -----------------------------------------------------------------
        % Print the nodal displacements
        function printResults(this)
            fprintf('******** NODAL DISPLACEMENTS ********\n');
            fprintf('\nNode       DX            DY\n');
            for nd = 1:this.nnodes
                fprintf('%2d     %10.3e    %10.3e\n',nd,...
                    this.U(this.ID(nd,1)),this.U(this.ID(nd,2)));
            end


            if ~isempty(this.enrDof)
                fprintf('\n******** ELEMENT ENRICHMENT *********\n');
                fprintf('\n\tFormulation:                 %s',  this.enhancementType);
                fprintf('\n\tEnrichment variable:         %s',  this.enrVar);
                fprintf('\n\tLevel of enrichment:         %s',  this.lvlEnrVar);
                fprintf('\n\tJump interpolation order:    %d',  this.jumpOrder);
                fprintf('\n\tConsider tangential stretch: %s',  mat2str(this.stretch(1)));
                fprintf('\n\tConsider normal stretch:     %s',  mat2str(this.stretch(2)));
                fprintf('\n\tApply static condensation:   %s\n',mat2str(this.staticCondensation));
                if (this.jumpOrder == 1) && strcmp(this.enrVar,'w')
                    fprintf('\nEl        wx1            wy1           wx2           wy2\n');
                elseif (this.jumpOrder == 1) && strcmp(this.enrVar,'alpha')
                    fprintf('\nEl        ax1            ay1           ax2           ay2\n');
                elseif (this.jumpOrder == 0) && strcmp(this.enrVar,'w')
                    fprintf('\nEl        wx1            wy1\n');
                elseif (this.jumpOrder == 0) && strcmp(this.enrVar,'alpha')
                    fprintf('\nEl        ax1            ay1\n');
                end
                for el = 1:this.nelem
                    if sum(this.IDenr(el,:)) > 0
                        if (this.jumpOrder == 0)
                            fprintf('%2d     %10.3e    %10.3e\n',nd,this.U(this.GLW{el}));
                        elseif (this.jumpOrder == 1)
                            fprintf('%2d     %10.3e    %10.3e    %10.3e    %10.3e\n',nd,this.U(this.GLW{el}));
                        end
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotMeshWithBC(this)
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh();
        end

        % -----------------------------------------------------------------
        % Plot the deformed mesh
        function plotDeformedMesh(this)

            this.updateResultVertices('Deformed');
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh();

        end

        %------------------------------------------------------------------
        % Update the result nodes coordinates of each element
        function updateResultVertices(this,configuration)
            
            for el = 1:this.nelem
                
                % Initialize the vertices array
                vertices = this.element(el).type.result.vertices0;

                % Get the updated vertices:
                if strcmp(configuration,'Deformed')

                    % Update the nodal displacement vector associated to the
                    % element. This displacement can contain the enhancement
                    % degrees of freedom.
                    this.element(el).type.ue = this.U(this.element(el).type.gle); 

                    % Update the vertices based on the displacement vector
                    % associated to the element
                    for i = 1:length(this.element(el).type.result.faces)
                        X = vertices(i,:);
                        u = this.element(el).type.displacementField(X);
                        vertices(i,:) = X + u';
                    end
                end
                this.element(el).type.result.setVertices(vertices);

                if isa(this.element(el).type,'EnrichedElement')
                    fractVertices = this.element(el).type.fracture.result.vertices0;
                    % Get the updated vertices:
                    if strcmp(configuration,'Deformed')
    
                        for i = 1:size(fractVertices,1)
                            X = fractVertices(i,:);
                            u = this.element(el).type.displacementField(X);
                            fractVertices(i,:) = X + u';
                        end
    
                    end
                    this.element(el).type.fracture.result.setVertices(fractVertices);
                end

            end
            
        end

        %------------------------------------------------------------------
        % Update the result nodes data of each element
        function updateResultVertexData(this,type)
            for el = 1:this.nelem
                % Update the nodal displacement vector associated to the
                % element. This displacement can contain the enhancement
                % degrees of freedom.
                this.element(el).type.ue = this.U(this.element(el).type.gle); 
                vertexData = zeros(length(this.element(el).type.result.faces),1);
                for i = 1:length(this.element(el).type.result.faces)
                    X = this.element(el).type.result.vertices(i,:);
                    if strcmp(type,'Model')
                        vertexData(i) = 0.0;
                    elseif strcmp(type,'Ux')
                        u = this.element(el).type.displacementField(X);
                        vertexData(i) = u(1);
                    elseif strcmp(type,'Uy')
                        u = this.element(el).type.displacementField(X);
                        vertexData(i) = u(2);
                    elseif strcmp(type,'Sx')
                        %
                    elseif strcmp(type,'Sy')
                        %
                    end
                end
                this.element(el).type.result.setVertexData(vertexData);
            end
        end

    end
end
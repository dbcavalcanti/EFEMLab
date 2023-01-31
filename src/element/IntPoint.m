%% IntPoint class
%
% This class defines an integration point object
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: January 2023
%%%
%
%% Class definition
classdef IntPoint < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        X               = [];  % Coordinates of the integration point in the natural coordinate system
        w               = 0.0; % Weight associated to the integration point
        strain          = [];  % Current strain vector
        stress          = [];  % Current stress vector
        statevar        = [];  % Current state variables vector
        strainOld       = [];  % Previous strain vector
        stressOld       = [];  % Previous stress vector  
        statevarOld     = [];  % Previous state variables vector
        constitutiveMdl = [];  % Constitutive model object
        anm             = '';  % Analysis model tag
        nVar            = 3;   % Dimension of the stress and strain vectors
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = IntPoint(X,w,anm,constitutiveMdl)
            if nargin > 0
                if strcmp(anm,'PlaneStress')
                    this.nVar = 3;
                elseif strcmp(anm,'PlaneStrain')
                    this.nVar = 3;
                elseif strcmp(anm,'Interface')
                    this.nVar = 2;
                end
                this.X               = X;
                this.w               = w;
                this.anm             = anm;
                this.strainOld       = zeros(this.nVar,1);
                this.stressOld       = zeros(this.nVar,1);
                this.statevarOld     = zeros(this.nVar,1);
                this.strain          = zeros(this.nVar,1);
                this.stress          = zeros(this.nVar,1);
                this.statevar        = zeros(this.nVar,1);
                this.constitutiveMdl = constitutiveMdl;
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        %  Compute the current strain vector
        function strainVct(this,dStrain)
            this.strain = this.strainOld + dStrain;
        end

        %------------------------------------------------------------------
        %  Update the current strain vector
        function updateStrainVct(this)
            this.strainOld = this.strain;
        end

        %------------------------------------------------------------------
        %  Compute the current stress vector
        function stressVct(this,dStrain)
            this.stress = this.constitutiveMdl.stressVct(dStrain,this);
        end

        %------------------------------------------------------------------
        %  Update the current stress vector
        function updateStressVct(this)
            this.stressOld = this.stress;
        end

        %------------------------------------------------------------------
        %  Get the current constitutive matrix
        function D = getConstitutiveMtrx(this,dStrain)
            D = this.constitutiveMdl.constitutiveMtrx(dStrain,this);
        end

        %------------------------------------------------------------------
        %  Get the current constitutive matrix
        function [stress,D] = constitutiveModel(this,dStrain)
            [stress,D] = this.constitutiveMdl.evalConstitutiveModel(dStrain,this);
        end

    end
end
%% Material_Elastic class
%
% This class defines an linear elastic stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef Material_Elastic < Material  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_Elastic(parameters, anm)
            this = this@Material('elastic',parameters, anm);
            this.nStVar = 0;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector
        function stress = stressVct(this, dStrain, pt)
            De      = this.constitutiveMtrx();
            stress  = De*(pt.strainOld + dStrain);
        end
        
        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix
        function De = constitutiveMtrx(this, ~, ~, ~, ~)

            % Elastic material properties
            E  = this.parameters(1);        % Young's modulus
            nu = this.parameters(2);        % Poisson coefficient

            if strcmp(this.anm,'PlaneStress')

                De = [ 1.0    nu    0.0;
                       nu    1.0    0.0;
                       0.0   0.0  (1-nu)/2.0 ] * E/(1.0 - (nu*nu));

            elseif strcmp(this.anm,'PlaneStrain')

                De = [ 1.0-nu    nu       0.0;
                         nu    1.0-nu     0.0;
                        0.0     0.0    (1-2.0*nu)/2.0 ] * E/(1.0 + nu)/(1.0 - 2.0*nu);

            else
                De = [];
            end
        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        function [stress,De] = evalConstitutiveModel(this,dStrain,pt)

            % Elastic material properties
            E  = this.parameters(1);        % Young's modulus
            nu = this.parameters(2);        % Poisson coefficient

            % Constitutive matrix
            if strcmp(this.anm,'PlaneStress')

                De = [ 1.0    nu    0.0;
                       nu    1.0    0.0;
                       0.0   0.0  (1-nu)/2.0 ] * E/(1.0 - (nu*nu));

            elseif strcmp(this.anm,'PlaneStrain')

                De = [ 1.0-nu    nu       0.0;
                         nu    1.0-nu     0.0;
                        0.0     0.0    (1-2.0*nu)/2.0 ] * E/(1.0 + nu)/(1.0 - 2.0*nu);

            else
                De = [];
            end

            % Stress vector
            stress  = De*(pt.strainOld + dStrain);

        end


    end
end
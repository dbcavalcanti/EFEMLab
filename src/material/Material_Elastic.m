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
        end
    end

    %% Public methods
    methods

        function stress = stressVct(this,dStrain)
            strain0 = this.point.strain0;
            De      = this.constitutiveMtrx();
            stress  = De*(strain0 + dStrain);
        end
        
        function De = constitutiveMtrx(this)
            if strcmp(this.anm,'PlaneStress')

            elseif strcmp(this.anm,'PlaneStrain')

            else
                De = [];
            end
        end

    end
end
%% MaterialInterface_Elastic class
%
% This class defines an elastic traction-displacement constitutive law
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
classdef MaterialInterface_Elastic < MaterialInterface      
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface_Elastic(parameters, penal)
            this = this@MaterialInterface('elastic', parameters, penal);
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the elastic traction vector
        %
        % Input:
        %   w0: current jump displacement vector in the local system [w0s, w0n]
        %   dw: increment of the jump displacement vector in the local
        %       system
        %
        % Output:
        %   t : traction vector in the local system [ts, tn]
        %
        function t = stressVct(this, dw, pt)

            Te = this.constitutiveMtrx(dw, pt);
            t  = Te*(pt.strainOld + dw);
            
        end

        %------------------------------------------------------------------
        % Compute the constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : elastic constitutive matrix
        %
        function Te = constitutiveMtrx(this, dw, pt)

            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

        end

        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : elastic constitutive matrix
        %
        function Te = elasticConstitutiveMtrx(this,w)

            % Penalization coefficient
            cs = 1.0;
            cn = 1.0;
            if this.penal == true
                if w(1) < 0.0
                    cs = 1.0e6;
                end
                if w(2) < 0.0
                    cs = 1.0e6;
                end
            end

            % Elastic material properties
            ks = this.parameters(1);
            kn = this.parameters(2);
            
            % Constitutive matrix
            Te = [ cs*ks   0.0;
                    0.0   cn*kn ];

        end

    end
end
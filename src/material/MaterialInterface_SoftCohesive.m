%% MaterialInterface_IsotropicDamage class
%
% This class defines an traction-displacement constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: February 2023
%%%
%
%% Class definition
classdef MaterialInterface_SoftCohesive < MaterialInterface_Elastic       
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface_SoftCohesive(parameters, penal)
            this = this@MaterialInterface_Elastic(parameters, penal);

            % This isotropic scalar damage model requires two state
            % variables, the maximum shear jump and the maximum (positive)
            % normal jump
            this.nStVar = 1; 

        end
 
    end

    %% Abstract methods
    methods(Abstract)

        % Compute the equivalent traction 
        equivalentTraction(this,pt,weq);

    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the traction vector
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

            % Secant constitutive matrix
            Tsec = this.secantConstitutiveMtrx(pt.strainOld + dw);

            % Stress vector
            t  = Tsec*(pt.strainOld + dw);
            
        end

        %------------------------------------------------------------------
        % Compute the tangent constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : tangent constitutive matrix
        %
        function T = constitutiveMtrx(this, dw, pt)

            % Elastic constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Scalar damage
            d = this.scalarDamage(dw, pt);

            % Tangent constitutive matrix
            T = (1.0 - d)*Te;

        end

        %------------------------------------------------------------------
        % Compute the secant constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : secant constitutive matrix
        %
        function Tsec = secantConstitutiveMtrx(this, dw, pt)

            % Elastic constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Scalar damage
            d  = this.scalarDamage(dw, pt);

            % Secant constitutive matrix
            Tsec = (1.0 - d)*Te;

        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        function [t,T] = evalConstitutiveModel(this,dw,pt)

            % Current total jump
            w = pt.strainOld + dw;

            % Compute the equivalent jump
            weq = this.evaluateStateVariable(w,pt);

            % Equivalent traction

            % Evaluate the damage criterion
            f = this.evaluateDamageCriterion(weq,pt);

            % Elastic constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Scalar damage
            [d,DdDweq]  = this.scalarDamage(w, pt);

            % Secant constitutive matrix
            Tsec = (1.0 - d)*Te;

            % Stress vector
            t  = Tsec * w;

            % Part of the tangent stiffness matrix associated with the
            % damage
            Td = this.damageConstitutiveMatrix(w,DdDweq);

            % Tangent constitutive matrix
            T = Tsec - Td;

        end

        %------------------------------------------------------------------
        % Update the scalar internal variables and the current maximum
        % values of the normal and shear jumps.
        function weq = evaluateStateVariable(this,w,pt)
            
            % Shear contribution factor
            beta = this.parameters(5);

            % Equivalent displacement jump
            weq = sqrt(w(2)*w(2) + beta*beta * w(1)*w(1));

            % Save the new state variables
            pt.statevar(1) = max(weq, pt.statevar(1));

        end

        %------------------------------------------------------------------
        % Evaluate the damage evolution criterion
        function f = evaluateDamageCriterion(~,weq,pt)

            % Evaluate the damage criterion
            f = weq - pt.statevar(1);

        end

        %------------------------------------------------------------------
        function Td = damageConstitutiveMatrix(this,w,DdDweq)

            % Get the material parameters
            ks = this.parameters(1); 
            kn = this.parameters(2);

            % Shear contribution factor
            beta = this.parameters(5);

            % Jump shear and normal components
            ws = w(1); wn = w(2);

            % Constitutive matrix
            Td = DdDweq * [ kn*wn  kn*wn*beta*sign(ws) ;
                            ks*ws  ks*ws*beta*sign(ws) ];
        end

    end

end
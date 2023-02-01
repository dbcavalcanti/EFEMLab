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
classdef MaterialInterface_IsotropicDamage < MaterialInterface_Elastic       
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface_IsotropicDamage(parameters, penal)
            this = this@MaterialInterface_Elastic(parameters, penal);

            % This isotropic scalar damage model requires two state
            % variables, the maximum shear jump and the maximum (positive)
            % normal jump
            this.nStVar = 2; 

        end
 
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
        % Damage scalar value. 
        % Considers a exponential evolution law
        function [d,DdDweq] = scalarDamage(this, w, pt)

            % Get the material parameters
            kn = this.parameters(2);
            ft = this.parameters(3);
            Gf = this.parameters(4);

            % Damage threshold
            w0 = ft / kn;

            % Compute the equivalent displacement jump
            weq0 = this.evaluateStateVariable([0.0,0.0],pt);
            weq  = this.evaluateStateVariable(w,pt);

            % Before damage initiated
            if (weq - w0) <= 0.0

                % Scalar damage
                d = 0.0;

                % Derivative of the scalar damage wrt the variable weq
                DdDweq = 0.0;

            % Damage initiated, but did not evolved
            elseif (weq - weq0) <= 0.0

                % Scalar damage
                d =  1.0 - w0/weq * exp(-ft*(weq - w0)/Gf); 

                % Derivative of the scalar damage wrt the variable weq
                DdDweq = 0.0;

            % Damage initiated and evolved
            else

                % Scalar damage
                d =  1.0 - w0/weq * exp(-ft*(weq - w0)/Gf);  
    
                % Derivative of the scalar damage wrt the variable weq
                DdDweq =  (ft*weq + Gf)* w0/(weq * weq * Gf) * exp(-ft*(weq - w0)/Gf); 

            end

        end

        %------------------------------------------------------------------
        % Update the scalar internal variables and the current maximum
        % values of the normal and shear jumps.
        function weq = evaluateStateVariable(this,w,pt)
            
            % Compute the maximum shear and normal jumps in the history of
            % this integration point
            wsMax = max(pt.statevarOld(1),abs(w(1)));
            wnMax = max(pt.statevarOld(2),max(w(2),0.0));

            % Shear contribution factor
            beta = this.parameters(5);

            % Equivalent displacement jump
            weq = wnMax + beta * wsMax;

            % Save the new state variables
            pt.statevar(1) = wsMax;
            pt.statevar(2) = wnMax;

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
%% Material_OrthotropicDamage class
%
% This class defines an orthotropic d+/d- damage model with memory constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef Material_OrthotropicDamage < Material  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_OrthotropicDamage(parameters, anm)
            this = this@Material('orthoDamage',parameters, anm);
            this.nStVar = 4;        % For plane strain problems
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector
        function stress = stressVct(this, dStrain, pt)

            % Secant constitutive matrix
            Dsec = this.secantConstitutiveMtrx(pt.strainOld + dStrain);

            % Stress vector
            stress  = Dsec*(pt.strainOld + dStrain);

        end
        
        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix (plane strain)
        function De = constitutiveMtrx(this, ~, ~, ~, ~)

            % Elastic material properties
            E  = this.parameters(1);        % Young's modulus
            nu = this.parameters(2);        % Poisson coefficient

            % Plane strain elastic constitutive matrix
            De = [ 1.0-nu    nu       0.0;
                     nu    1.0-nu     0.0;
                    0.0      0.0    (1-2.0*nu)/2.0 ] * E/(1.0 + nu)/(1.0 - 2.0*nu);
        end

        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix (plane strain)
        function De = elasticConstitutiveMtrx(this)

            % Elastic material properties
            E  = this.parameters(1);        % Young's modulus
            nu = this.parameters(2);        % Poisson coefficient

            % Plane strain elastic constitutive matrix
            De = [ 1.0-nu    nu       0.0;
                     nu    1.0-nu     0.0;
                    0.0      0.0    (1-2.0*nu)/2.0 ] * E/(1.0 + nu)/(1.0 - 2.0*nu);
        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        function [stress,D] = evalConstitutiveModel(this,dStrain,pt)

            % Current total strain
            strain = pt.strainOld + dStrain;

            % Secant constitutive matrix
            Dsec = this.secantConstitutiveMtrx(strain, pt);

            % Stress vector
            stress = Dsec * strain;

            % Tangent constitutive matrix
            D = this.tangentConstitutiveMtrx();

        end

        %------------------------------------------------------------------
        % Compute the secant constitutive matrix
        function Dsec = secantConstitutiveMtrx(this,strain,pt)

            % Elastic constitutive matrix
            De = this.elasticConstitutiveMtrx();

            % Compute the matrix with the damage parameters
            Dam = this.damageParamMtrx(strain,pt);

            % Secant constitutive matrix
            Dsec = Dam .* De;

        end

        %------------------------------------------------------------------
        % Compute the matrix with the damage parameters associated wih each
        % component of the constitutive matrix
        function Dam = damageParamMtrx(strain,pt)

            % Compute the vector damage parameter associated to each
            % orthogonal direction
            d = this.damageVct(strain,pt);

            % Initialize the matrix Dam
            Dam = zeros(3,3);            % Particular to 2D problems

            % Fill the matrix Dam
            Dam(1,1) = 1.0 - d(1);
            Dam(1,2) = min([1.0 - d(1); 1.0 - d(2)]);
            Dam(2,1) = Dam(1,2);
            Dam(2,2) = 1.0 - d(2);
            Dam(3,3) = Dam(1,2);

        end

        %------------------------------------------------------------------
        % Compute the vector damage parameter associated to each orthogonal 
        % direction
        function d = damageVct(this,strain,pt)

            % Compute the effective stress vector
            effStress = this.effectiveStressVct(strain);

            % Rotate the effective stress vector
            effStress = this.rotateStressVct(effStress,pt);

            % Compute the equivalent tensile/compressive effective stress
            % vector
            eqStress = this.equivalentStressVct(effStress);

            % Update the tensile/compressive damage threshold
            this.evaluateStateVariable(eqStress,pt);

            % Update the damage vector
            d = this.damageLaw(pt,sign(effStress));

        end

        %------------------------------------------------------------------
        % Compute the effective stress vector
        function effStress = effectiveStressVct(this,strain)

            % Elastic constitutive matrix
            De = this.elasticConstitutiveMtrx();

            % Compute the effective stress vector
            effStress = De * strain;

        end

        %------------------------------------------------------------------
        % Rotate the stress vector based on a given orientation vector
        function stress = rotateStressVct(this,stress,pt)

            % Compute the orientation vector 
            theta = this.orientationVct(stress,pt);

            % Get the orientation vector components
            cs = theta(1); sn = theta(2);

            % Compute the rotation matrix
            R = [cs*cs   sn*sn   -2.0*cs*sn;
                 sn*sn   cs*cs    2.0*cs*sn;
                 cs*sn  -cs*sn   cs*cs-sn*sn];

            % Rotate the stress vector
            stress = R * stress;

        end

        %------------------------------------------------------------------
        % Compute the orientation vector 
        function theta = orientationVct(this,stress,pt)

            if pt.statevar(end) > 0.0
                % Rotate to the principal stresses direction
                theta = this.principalStressOrient(stress);
            else
                % Rotate to the material local system
                theta = pt.statevar(3:4);  % Specific for 2D problems
            end

        end

        %------------------------------------------------------------------
        % Compute the principal stress orientation vector
        function thetap = principalStressOrient(~,stress)

            if (abs(stress(1) - stress(2)) > 0.0)
                thetap = 0.5 * atan2(stress(3),stress(1) - stress(2));
            elseif (txy > 0.0)
                thetap = pi / 4.0;
            elseif (txy < 0.0)
                thetap = -pi / 4.0;
            else
                thetap = 0.0;
            end
            
            if( thetap < 0.0 )
                thetap = pi + thetap;
            end

        end

        %------------------------------------------------------------------
        % Compute the equivalent stress vector based on the Rankine
        % criterion.
        % This vector is assembled such that the first to components are
        % associated to the tension in each direction and the last one are
        % associated to the compression.
        function eqStress = equivalentStressVct(~,stress)

            % Initialize the equivalent stress vector (2D problems)
            eqStress = zeros(4,1);

            % Fill the equivalent stress vector based on the Rankine
            % criterion
            if stress(1) > 0.0
                eqStress(1) =  stress(1);
            else
                eqStress(3) = -stress(1);
            end
            if stress(2) > 0.0
                eqStress(2) =  stress(2);
            else
                eqStress(4) = -stress(2);
            end

        end

        %------------------------------------------------------------------
        % Update the damage threshold state variables at the integration
        % point
        function evaluateStateVariable(~,eqStress,pt)

            pt.statevar(1) = max(pt.statevar(1), eqStress(1));
            pt.statevar(2) = max(pt.statevar(2), eqStress(2));
            pt.statevar(3) = max(pt.statevar(3), eqStress(3));
            pt.statevar(4) = max(pt.statevar(4), eqStress(4));

        end

        %------------------------------------------------------------------
        % Evaluate the parabolic-exponential softening damage law
        function d = damageLaw(this,pt,signStress)

            % Get the material parameters
            E   = this.parameters(1);       % Young's modulus
            ft  = this.paramaters(3);       % Tensile strength
            fc  = this.parameters(4);       % Compressive strength
            Gft = this.parameters(5);       % Tensile fracture energy
            Gfc = this.parameters(6);       % Compressive fracture energy
            g0t = this.parameters(7);       % Tensile gamma_0
            g0c = this.parameters(8);       % Compressive gamma_0
            gpt = this.parameters(9);       % Tensile gamma_p
            gpc = this.parameters(10);      % Compressive gamma_p
            b   = this.parameters(11);      % Crack width

            % Derivated auxiliary material parameters
            f0t  = g0t * ft;
            f0c  = g0c * fc;
            fpt  = gpt * ft;
            fpc  = gpc * fc;
            Adt  = (fpt - ft)/ft;
            Adc  = (fpc - fc)/fc;
            Adbt = Adt * (fpt^3 - 3.0*fpt*f0t^2 + 2.0*f0t^3)/(3.0*ft*(fpt - f0t)^2);
            Adbc = Adc * (fpc^3 - 3.0*fpc*f0c^2 + 2.0*f0c^3)/(3.0*fc*(fpc - f0c)^2);

            % Softening modulus
            Hdt = 1.0 / ((2.0 * E * Gft)/(ft * ft * b) - fpt/ft - Adbt);
            Hdc = 1.0 / ((2.0 * E * Gfc)/(fc * fc * b) - fpc/fc - Adbc);

            % Get the threshold parameter
            r = pt.statevar(1:4);
            
            % Initialize the damage vector
            d = zeros(2,1);

            % Parabolic-exponential softening law
            if signStress(1) > 0.0
                if (r(1) >= f0t) && (r(1) < fpt)
                    d(1) = Adt * (ft/r(1)) * ((r(1) - f0t)/(fpt - f0t));
                elseif r(1) >= fpt
                    d(1) = 1.0 - (ft/r(1)) * exp(2.0 * Hdt * ((ft - r(1))/ft));
                end
            elseif signStress(1) < 0.0
                if (r(3) > f0c) && (r(3) < fpc)
                    d(1) = Adc * (fc/r(3)) * ((r(3) - f0c)/(fpc - f0c));
                elseif r(3) > fpc
                    d(1) = 1.0 - (fc/r(3)) * exp(2.0 * Hdc * ((fc - r(3))/fc));
                end
            end
            if signStress(2) > 0.0
                if (r(2) > f0t) && (r(2) < fpt)
                    d(2) = Adt * (ft/r(2)) * ((r(2) - f0t)/(fpt - f0t));
                elseif r(2) > fpt
                    d(2) = 1.0 - (ft/r(2)) * exp(2.0 * Hdt * ((ft - r(2))/ft));
                end
            elseif signStress(2) < 0.0
                if (r(4) > f0c) && (r(4) < fpc)
                    d(2) = Adc * (fc/r(4)) * ((r(4) - f0c)/(fpc - f0c));
                elseif r(4) > fpc
                    d(2) = 1.0 - (fc/r(4)) * exp(2.0 * Hdc * ((fc - r(4))/fc));
                end
            end
            

        end

    end
end
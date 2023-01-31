%% TractionConstitutiveLaw_IsotropicDamage class
%
% This class defines an traction-displacement constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
%
%% Class definition
classdef TractionConstitutiveLaw_IsotropicDamage < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        k        = 0.0;      % Stiffness (isotropic material)
        d        = 0.0;      % Scalar damage
        kappa    = 0.0;      % Internal variable (largest jump history)
        kappa0   = 0.0;      % Internal variable relative initial value
        DdDkappa = 0.0;      % Derivative of the scalar damage wrt to kappa
        wn_max   = 0.0;      % Maximum normal jump (only positive values)
        ws_max   = 0.0;      % Maximum shear jump (absolute values)
        Gf       = 0.0;      % Fracture energy
        ft0      = 0.0;      % Initial tensile strength
        beta     = 0.0;      % Contribution of the shear jump component
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = TractionConstitutiveLaw_IsotropicDamage(k,d)
            if (nargin > 0)
                this.k = k;
                this.d = d;
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Update the scalar internal variables and the current maximum
        % values of the normal and shear jumps.
        function updateKappa(this,wn,ws)
            this.wn_max = max(this.wn_max,max(wn,0.0));
            this.ws_max = max(this.ws_max,abs(ws));
            this.kappa  = this.wn_max + this.beta * this.ws_max;
        end

        %------------------------------------------------------------------
        % Damage scalar value. 
        % Considers a exponential evolution law
        function updateScalarDamage(this)
            this.d =  1 - this.kappa0/this.kappa * exp(-this.ft0*(this.kappa - this.kappa0)/this.Gf);       
            this.DdDkappa =  (this.ft0*this.kappa + this.Gf)* this.kappa0/(this.kappa * this.kappa * this.Gf) * exp(-this.ft0*(this.kappa - this.kappa0)/this.Gf); 
        end

        %------------------------------------------------------------------
        % Update the scalar internal variables and the current maximum
        % values of the normal and shear jumps.
        function evaluateLoadingFnc(this,wn,ws)
            this.f  = max(wn,0.0) + this.beta * abs(ws) - this.kappa;
        end

        %------------------------------------------------------------------
        % Computes the elastic constitutive matrix
        function Tel = elasticConstitutiveMtrx(this)
            Tel = this.k*eye(2);
        end

        %------------------------------------------------------------------
        % Computes the elastic constitutive matrix
        function Td = damagedConstitutiveMtrx(this,wn,ws)
            Td  = zeros(2,2);
            Td(1,1) = -this.k*(this.d + this.DdDkappa*wn);
            Td(1,2) = -this.k * this.DdDkappa * this.beta * wn * sign(ws);
            Td(2,1) = -this.k * this.DdDkappa * ws;
            Td(2,2) = -this.k*(this.d + this.DdDkappa*ws*sign(ws));
        end

        %------------------------------------------------------------------
        % Computes the tangent constitutive matrix
        function tangentConstitutiveMtrx(this,wn,ws)

            % Update the internal variables
            this.updateKappa(wn,ws);
            this.updateScalarDamage();

            % Compute the elastic stiffness matrix
            Tel = elasticConstitutiveMtrx(this);

            % Compute the damaged part of the stiffness matrix
            Td = damagedConstitutiveMtrx(wn,ws);
            
            % Compute the tangent stiffness matrix
            this.T = Tel + Td;
        end


    end

end
%% Material class
%
% This class defines an abstract stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef Material < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        model       = 'elastic';   
        parameters  = [];
        anm         = 'PlaneStress';
        point       = [];
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material(model,parameters, point)
            this.model      = model;
            this.parameters = parameters;
            this.anm        = anm;
            this.point      = point;
        end
    end

    %% Abstract methods
    methods(Abstract)

        % Compute the stress vector
        stress = stressVct(this);

        % Compute the constitutive matrix
        De = constitutiveMtrx(this);
        
    end
end
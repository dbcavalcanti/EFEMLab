%% MaterialInterface class
%
% This class defines an abstract traction-displacement constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: January 2023
%%%
%% Class definition
classdef MaterialInterface < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        model       = 'elastic';   
        parameters  = [];
        penal       = false;   
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface(model, parameters, penal)
            if (nargin > 0)
                this.model      = model;
                this.parameters = parameters;
                this.penal      = penal;
            end
        end
    end
    
    %% Public methods
    methods(Abstract)

        %------------------------------------------------------------------
        % Compute the stress vector
        stress = tractionVct(this);

        %------------------------------------------------------------------
        % Compute the constitutive matrix
        De = constitutiveMtrx(this);

    end
end
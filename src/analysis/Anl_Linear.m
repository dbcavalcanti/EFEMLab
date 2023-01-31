%% Anl_Nonlinear Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with linear-elastic analysis.
%
classdef Anl_Linear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Linear()
            this = this@Anl('Linear');
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function process(this,mdl)

            % Assemble elastic stiffness matrix and check model stability
            K = mdl.globalStiffnessMtrx(mdl.U); 
            
            % Solve linear-elastic analysis
            [mdl.U, mdl.F] = this.solveSystem(mdl,K,mdl.F,mdl.U);

        end
    end
    
end
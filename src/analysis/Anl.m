%% Anl class
%
% This in an abstract class that defines the analysis object
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
classdef Anl < handle

    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        type     = [];          
    end

    %% Constructor method
    methods
        function this = Anl(type)
            if nargin > 0
                this.type = type;
            end
        end
    end
    methods
        %------------------------------------------------------------------
        % Partition and solve a linear system of equations.
        %  f --> free d.o.f. (natural B.C. - unknown) 
        %  c --> constrained d.o.f. (essential B.C. - known) 
        %
        % [ Kff Kfs ] * [ Uf ] = [ Fext ]
        % [ Ksf Kss ]   [ Us ] = [   R  ]
        %
        function [U,Fext] = solveSystem(~,mdl,K,Fext,U)

            if nargin < 5
                U = zeros(mdl.nTotDofs,1);
            end

            % Partition system of equations
            freedof  = [1:mdl.ndoffree,mdl.enrFreeDof'];
            fixeddof = (1+mdl.ndoffree):mdl.ndof;
            Kff      = K(freedof, freedof);
            Kfs      = K(freedof, fixeddof);
            Ksf      = K(fixeddof,freedof);
            Kss      = K(fixeddof,fixeddof);
            Ff       = Fext(freedof);
            Fs       = Fext(fixeddof);
            Us       = U(fixeddof); 
            
            % Solve the system of equilibrium equations
            Uf = (Kff) \ (Ff - Kfs * Us);
            
            % Compute the reaction forces
            Fs = -Fs + Ksf * Uf + Kss * Us;
            Fext(fixeddof) = Fs;
            
            % Displacement vector
            U(freedof)  = Uf;
            
        end

    end
end
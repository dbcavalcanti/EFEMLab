%% Anl_Nonlinear Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with nonlinear analysis.
%
classdef Anl_Nonlinear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        method     = 0;   % flag for solution method
        incr_type  = 0;   % flag for type of increment size adjustment
        iter_type  = 0;   % flag for type of iteration strategy
        increment  = 0;   % initial increment of load ratio
        max_lratio = 0;   % limit value of load ratio
        max_step   = 0;   % maximum number of steps
        max_iter   = 0;   % maximum number of iterations in each step
        trg_iter   = 0;   % desired number of iterations in each step
        tol        = 0;   % numerical tolerance for convergence
        ctrlNode   = 0;   % control node (for displacement control method)
        ctrlDof    = 0;   % control dof (for displacement control method)
        incrSign   = 0;   % sign of the increment of displacement
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl_Nonlinear()
            anl = anl@Anl('Nonlinear');
            
            % Default analysis options
            anl.method     = c.LOAD_CONTROL;
            anl.incr_type  = c.CONSTANT;
            anl.iter_type  = c.STANDARD;
            anl.increment  = 0.01;
            anl.max_lratio = 1.0;
            anl.max_step   = 1000;
            anl.max_iter   = 100;
            anl.trg_iter   = 3;
            anl.tol        = 0.000001;
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function status = process(anl,sim)
            c = Constants();
            status = 1;
            mdl = sim.mdl;
            res = sim.res;
            tic
            
            % Check for any free d.o.f
            if (mdl.neqf == 0)
			    status = 0;
                fprintf('Status: Model with no free degree-of-freedom!\n');
                return;
            end
            
            % Initialize results
            res.steps = 0;
            res.lbd   = zeros(anl.max_step+1,1);
            res.U     = zeros(mdl.neq,anl.max_step+1);
            res.niter = 0;
            
            % Initialize data for first step
            step  = 0;  % step number
            lbd   = 0;  % total load ratio (lambda)
            sign  = 1;  % sign of predicted increment of load ratio
            
            % Initialize vector of total nodal displacements and reference load vector
            U    = zeros(mdl.neq,1);
            Pref = zeros(mdl.neq,1);

            % Initialize vector of total increment displacement
            D_U = zeros(mdl.neq,1);
            
            % Add contributions of nodal forces, element loads, and presc. displ. to reference load vector
            Pref = mdl.addNodalLoad(Pref);
            Pref = mdl.addEquivLoad(Pref);
            Pref = mdl.addPrescDispl(Pref);
            
            %==========================================================================
            % Start incremental process
            while (step < anl.max_step)
                step = step + 1;
                
                % Tangent stiffness matrix
                Kt = mdl.gblTangStiffMtx(anl.tang_mtx);
                
                % Tangent increment of displacements for predicted solution
                d_Up0 = anl.solveSystem(mdl,Kt,Pref,true,false);
                
                if (step == 1)
                    % Initial increment of load ratio for predicted solution
                    if (anl.method == c.DISPL_CONTROL)
                        d_lbd0 = anl.predictedIncrement(anl,mdl,sign,1,1,0.0,0.0,D_U,d_Up0,Pref,c);
                    else
                        d_lbd0 = anl.increment;
                    end
                    
                    % Set previous tangent increment of displacements as current increment
                    d_Up0_old = d_Up0;
                    
                    % Store squared value of the norm of tangent increment of displacements
                    n2 = d_Up0(1:mdl.neqf)'*d_Up0(1:mdl.neqf);
                else
                    % Generalized Stiffness Parameter
                    GSP = n2/(d_Up0(1:mdl.neqf)'*d_Up0_old(1:mdl.neqf));
                    
                    % Adjust increment sign
                    if (GSP < 0)
                        sign = -sign;
                    end
                    
                    % Adjustment factor of increment size
                    if (anl.incr_type == c.CONSTANT)
                        J = 1;
                    elseif (anl.incr_type == c.ADJUSTED)
                        J = sqrt(anl.trg_iter/iter);
                    end
                    
                    % Predicted increment of load ratio
                    d_lbd0 = anl.predictedIncrement(anl,mdl,sign,J,GSP,D_lbd,d_lbd0,D_U,d_Up0,Pref,c);
                end
                
                % Limit increment of load ratio to make total load ratio smaller than maximum value
                if ((anl.max_lratio > 0.0 && lbd + d_lbd0 > anl.max_lratio) ||...
                    (anl.max_lratio < 0.0 && lbd + d_lbd0 < anl.max_lratio))
                    d_lbd0 = anl.max_lratio - lbd;
                end
                
                % Increments of load ratio and displacements for predicted solution
                d_lbd = d_lbd0;
                d_U0  = d_lbd0 * d_Up0;
                d_U   = d_U0;
                
                % Initialize incremental values of load ratio and displacements for current step
                D_lbd = d_lbd;
                D_U   = d_U;
                
                % Update total values of load ratio and displacements
                lbd = lbd + d_lbd;
                U   = U   + d_U;
                
                %----------------------------------------------------------------------
                % Start iterative process
                iter = 1;
                conv = 0;
                while (conv == 0 && iter <= anl.max_iter)
                    
                    % Reset reference load vector
                    Pref = mdl.F;
                    
                    % Vector of external and internal forces
                    Fext = lbd * Pref;
                    Fint = mdl.globalInternalForceVct(dU);
                    
                    % Vector of unbalanced forces
                    R = Fext - Fint;
                    
                    % Check convergence
                    unbNorm = norm(R(1:mdl.neqf));
                    forNorm = norm(Pref(1:mdl.neqf));
                    conv = (unbNorm == 0 || forNorm == 0 || unbNorm/forNorm < anl.tol);
                    if conv == 1
                        break;
                    end
                    
                    % Start/keep corrector phase
                    iter = iter + 1;
                    
                    % Tangent stiffness matrix
                    K = mdl.globalStiffnessMtrx(dU); 
                    
                    % Tangent and residual increments of displacements
                    d_Up = anl.solveSystem(mdl,K,Pref);
                    d_Ur = anl.solveSystem(mdl,K,R);
                    
                    % Corrected increment of load ratio
                    d_lbd = anl.correctedIncrement(anl,mdl,d_lbd0,D_lbd,d_Up0_old,d_U0,d_Up,d_Ur,D_U,Pref,R,c);
                    if (~isreal(d_lbd))
                        conv = -1;
                        break;
                    end
                    
                    % Corrected increment of displacements
                    d_U = d_lbd * d_Up + d_Ur;
                    
                    % Increments of load ratio and displacements for current step
                    D_lbd = D_lbd + d_lbd;
                    D_U   = D_U   + d_U;
                    
                    % Total values of load ratio and displacements
                    lbd = lbd + d_lbd;
                    U   = U   + d_U;
                end
                %----------------------------------------------------------------------
                % Check for convergence fail or complex value of increment
                if (conv == 0)
                    status = (step > 1);
                    fprintf('Status: Convergence not achieved!\n');
                    return;
                elseif (conv == -1)
                    status = (step > 1);
                    fprintf('Status: Unable to compute load increment!\n');
                    return;
                end
                
                % Store step results
                res.steps = step;
                res.lbd(step+1) = lbd;
                res.U(:,step+1) = U;
                res.niter = res.niter + iter;
                
                % Store predicted tangent increment of displacements for next step
                if (step ~= 1)
                    d_Up0_old = d_Up0;
                end
                
                % Check if maximum load ratio was reached
                if ((anl.max_lratio >= 0 && lbd >= 0.999*anl.max_lratio) ||...
                    (anl.max_lratio <= 0 && lbd <= 0.999*anl.max_lratio))
                    break;
                end
            end
            %==========================================================================
            
            % Clean unused steps
            if (step < anl.max_step)
                res.lbd = res.lbd(1:step+1);
                res.U   = res.U(:,1:step+1);
            end

        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Compute inrement of load ratio for the predicted solution
        % (first iteration)
        function d_lbd0 = predictedIncrement(anl,mdl,sign,J,GSP,D_lbd,d_lbd0,D_U,d_Up0,Pref,c)
            % Extract free d.o.f. components
            Pref  = Pref(1:mdl.neqf);
            D_U   = D_U(1:mdl.neqf);
            d_Up0 = d_Up0(1:mdl.neqf);
            
            % EULER: GSP criteria
            if (anl.method == c.EULER)
                if (anl.incr_type == c.CONSTANT)
                    d_lbd0 = anl.increment;
                elseif (anl.incr_type == c.ADJUSTED)
                    d_lbd0 = sqrt(abs(GSP)) * anl.increment;
                end
                
            % LCM: Load Increment
            elseif (anl.method == c.LOAD_CONTROL)
                d_lbd0 = J * abs(d_lbd0);

            % DCM: Displacement Increment
            elseif (anl.method == c.DISPL_CONTROL)
                dof = mdl.ID(mdl.anm.gla==anl.ctrlDof,anl.ctrlNode);
                d_lbd0 = J * anl.incrSign * anl.increment / d_Up0(dof);
                return  % The change of the sign must not be applied
                
            % WCM: Work Increment
            elseif (anl.method == c.WORK_CONTROL)
                d_lbd0 = J * sqrt(abs((D_lbd*Pref'*D_U)/(Pref'*d_Up0)));
                
            % ALCM_FNP: Cylindrical Arc-Length Increment
            elseif (anl.method == c.ARC_LENGTH_FNP)
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ALCM_UNP: Cylindrical Arc-Length Increment
            elseif (anl.method == c.ARC_LENGTH_UNP)
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ALCM_CYL: Cylindrical Arc-Length Increment
            elseif (anl.method == c.ARC_LENGTH_CYL)
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ALCM_SPH: Spherical Arc-Length Increment
            elseif (anl.method == c.ARC_LENGTH_SHP)
                d_lbd0 = J * sqrt((D_U'*D_U + D_lbd^2*(Pref'*Pref)) / (d_Up0'*d_Up0 + Pref'*Pref));
                
            % MNCM: Cylindrical Arc-Length Increment
            elseif (anl.method == c.MINIMUM_NORM)
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % ORCM: Cylindrical Arc-Length Increment
            elseif (anl.method == c.ORTHOGONAL_RES)
                d_lbd0 = J * sqrt((D_U'*D_U)/(d_Up0'*d_Up0));
                
            % GDCM: GSP criteria
            elseif (anl.method == c.GENERAL_DISPL)
                d_lbd0 = J * sqrt(abs(GSP)) * anl.increment;
            end
            
            % Apply increment sign
            d_lbd0 = sign * d_lbd0;
        end
        
        %--------------------------------------------------------------------------
        % Compute inrement of load ratio for the corrected solutions
        % (iterations to correct predicted solution).
        function d_lbd = correctedIncrement(anl,mdl,d_lbd0,D_lbd,d_Up0,d_U0,d_Up,d_Ur,D_U,Pref,R,c)
            % Extract free d.o.f. components
            d_Up0 = d_Up0(1:mdl.neqf);
            d_U0  = d_U0(1:mdl.neqf);
            d_Up  = d_Up(1:mdl.neqf);
            d_Ur  = d_Ur(1:mdl.neqf);
            D_U   = D_U(1:mdl.neqf);
            Pref  = Pref(1:mdl.neqf);
            R     = R(1:mdl.neqf);
            
            % LCM
            if (anl.method == c.LOAD_CONTROL)
                d_lbd = 0;
            
            % DCM
            elseif (anl.method == c.DISPL_CONTROL)
                dof = mdl.ID(mdl.anm.gla==anl.ctrlDof,anl.ctrlNode);
                d_lbd = -d_Ur(dof)/d_Up(dof);
                
            % WCM
            elseif (anl.method == c.WORK_CONTROL)
                d_lbd = -(Pref'*d_Ur)/(Pref'*d_Up);
                
            % ALCM_FNP
            elseif (anl.method == c.ARC_LENGTH_FNP)
                d_lbd = -(d_Ur'*d_U0)/(d_Up'*d_U0 + d_lbd0*(Pref'*Pref));
                
            % ALCM_UNP
            elseif (anl.method == c.ARC_LENGTH_UNP)
                d_lbd = -(d_Ur'*D_U)/(d_Up'*D_U + D_lbd*(Pref'*Pref));
                
            % ALCM_CYL
            elseif (anl.method == c.ARC_LENGTH_CYL)
                a = d_Up'*d_Up;
                b = d_Up'*(d_Ur + D_U);
                c = d_Ur'*(d_Ur + 2*D_U);
                s = sign(D_U'*d_Up);
                
                d_lbd = -b/a + s*sqrt((b/a)^2 - c/a);
                
            % ALCM_SPH
            elseif (anl.method == c.ARC_LENGTH_SHP)
                a = d_Up'*d_Up + Pref'*Pref;
                b = d_Up'*(d_Ur + D_U) + D_lbd*(Pref'*Pref);
                c = d_Ur'*(d_Ur + 2*D_U);
                s = sign(D_U'*d_Up);
                
                d_lbd = -b/a + s*sqrt((b/a)^2 - c/a);
                
            % MNCM
            elseif (anl.method == c.MINIMUM_NORM)
                d_lbd = -(d_Up'*d_Ur)/(d_Up'*d_Up);
                
            % ORCM
            elseif (anl.method == c.ORTHOGONAL_RES)
                d_lbd = -(R'*D_U)/(Pref'*D_U);
                
            % GDCM
            elseif (anl.method == c.GENERAL_DISPL)
                d_lbd = -(d_Up0'*d_Ur)/(d_Up0'*d_Up);
            end
        end
        
    end
end
classdef AvancoTemporal

    properties (Constant)
        tempo_final = 10;
        dt = 0.01; 
    end
    
    properties (Constant, Access = public)
        num_steps = ceil(AvancoTemporal.tempo_final / AvancoTemporal.dt);
        tempo_final_real = AvancoTemporal.num_steps * AvancoTemporal.dt;
    end

    methods (Static)
         
        % function [rhs_u, rhs_v, rhs_k, rhs_epsilon] = RK4Step(u, v, k, epsilon)
        %     [rhs_u, rhs_v] = NavierStokesRHS.CalculateRHS(u, v, k, epsilon);
        %     [rhs_k] = kineticEnergyRHS.CalculateRHS(u, v, k, epsilon);
        %     [rhs_epsilon] = epsilonRHS.CalculateRHS(u, v, k, epsilon);
        % end
        
    end

end

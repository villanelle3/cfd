classdef Canal
    properties (Constant)
        % Parâmetros do canal
        Lx = 10;                % Comprimento do canal
        Ly = 1;                 % Largura do canal
        nx = 8;                 % Número de pontos na direção do comprimento
        ny = 16;                % Número de pontos na direção da parede
        Re = 40000;             % Número de Reynolds
        nu = 1/Canal.Re;        % Viscosidade cinemática
        uinlet = 1.0;

        % Constantes do modelo
        rho = 1; 
        c_mu = 0.09;
        sigma_k = 1;
        sigma_eps = 1.22;
        c_1eps = 1.44;
        c_2eps = 1.92;
        omegaeps = 0.3;
        omegak = 0.3;   
        kappa = 0.41;   
        B = 5.5;
        ti = 0.1;               % Intensidade turbulenta
    end
    
    properties (Constant, Access = public)
        x = linspace(0, Canal.Lx, Canal.nx);   
        y = linspace(0, Canal.Ly, Canal.ny);
        dx = Canal.Lx / (Canal.nx - 1);         % Espaçamento de grade na direção do comprimento (x)
        dy = Canal.Ly / (Canal.ny - 1);         % Espaçamento de grade na direção da largura (y)
    end
    
    methods (Static)
        function [u, v, nu_til, mu_t] = AplicarCondicaoInicial()
            % Método para aplicar as condições iniciais
            
            % Inicializar matrizes para as variáveis do escoamento
            u = Canal.uinlet * ones(Canal.nx, Canal.ny); % Componente da velocidade na direção do escoamento 
            v = zeros(Canal.nx, Canal.ny);               % Componente da velocidade na direção da largura
            
            % Garantir que as primeiras e últimas linhas de u e v representem as paredes
            u(:, 1) = 0;        % Velocidade na parede inferior
            u(:, end) = 0;      % Velocidade na parede superior
            v(:, 1) = 0;        % Velocidade na parede inferior na direção da largura
            v(:, end) = 0;      % Velocidade na parede superior na direção da largura
            
            % Calcular energia cinética turbulenta (k)
            k = 2/3 * (u .* Canal.ti).^2;
            
            % Calcular taxa de dissipação turbulenta (epsilon)
            Ls = 0.07 * Canal.Ly;
            epsilon = Canal.c_mu^(3/4) * k.^(3/2) ./ Ls;
            
            % Calcular viscosidade turbulenta (mu_t)
            mu_t = Canal.calcularMut(k, epsilon);

            nu_til = 0.01 * Canal.nu * ones(Canal.nx, Canal.ny);
            nu_til(:, 1) = 0;        % Velocidade na parede inferior
            nu_til(:, end) = 0;      % Velocidade na parede superior
        end

        function mu_t = calcularMut(k, epsilon)
            % Calcular viscosidade turbulenta (mu_t)
            mu_t = Canal.rho * Canal.c_mu .* k.^2 ./ epsilon;
            
            % Substituir valores NaN e Inf em mu_t por 0
            mu_t(isnan(mu_t) | isinf(mu_t)) = 0;
        end

        function [kx, ky, k2] = CalcularNumerosOnda()
            nx = Canal.nx;
            ny = Canal.ny;
            dx = Canal.dx;
            dy = Canal.dy;
            % Pré-alocação dos vetores de número de onda na direção "x" e "y"
            kx = zeros(1, nx);
            ky = zeros(1, ny);
            k2 = zeros(nx, ny);
        
            % Definição dos vetores de número de onda na direção "x"
            for i = 1:nx/2+1
                kx(i) = 2.0*pi*(i-1)/(nx*dx);
            end
        
            for i = nx/2+2:nx
                kx(i) = 2.0*pi*(i-1-nx)/(nx*dx);
            end
        
            % Definição dos vetores de número de onda na direção "y"
            for j = 1:ny/2+1
                ky(j) = 2.0*pi*(j-1)/(ny*dy);
            end
        
            for j = ny/2+2:ny
                ky(j) = 2.0*pi*(j-1-ny)/(ny*dy);
            end
        
            % Constante pequena para estabilidade numérica
            small = 1.0e-30;
        
            % Cálculo do quadrado do número de onda k2
            for j = 1:ny
                for i = 1:nx
                    k2(i,j) = kx(i)^2 + ky(j)^2 + small;
                end
            end
        end
    end
end

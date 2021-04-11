function T_new = SIMO_compute_Tnew(T_old, z, dt, TURB_DIF, CONV_MIX, LAT, u_star, z0, Q, Sn_0)

    %% SIMO v1.0: Simplified model of the vertical temperature profile in a small warm monomictic lake
    %% Last update 11.04.2021 
    %% Author: Kristina Šaroviæ
    %% 
    % A simple 1-D energy budget model for prediction of vertical temperature 
    % profile in small, monomictic lake, that is forced by reduced number of 
    % input meteorological variables is proposed. The model estimates the net
    % heat flux and thermal diffusion using only routinely measured hourly mean
    % meteorological variables (namely, the air temperature, relative humidity, 
    % atmospheric pressure, wind speed, and precipitation), hourly mean 
    % ultraviolet B radiation (UVB), and climatological monthly mean cloudiness data.    
    %
    % This function performs the integration step to calculate the temperature 
    % in the next time step.
    %
    %% INPUT %% -----------------------------------------------------------
        % T_old       % Water tempearture profile before the integation [ C ]
        % z           % Depths of the integration points                [ m ]
        % dt          % Time step (dt=1h=3600s recomended)              [ s ]     
        % TURB_DIF    % Turbulent difusion flag (included/neglected)    [ 1/0 ]
        % CONV_MIX    % Convective mixing  flag (included/neglected)    [ 1/0 ]
        % LAT         % Location latitude                               [ LAT ]
        % u_star      % Air shear velocity                              [ m s-1 ]
        % z0          % Roughness length                                [ m ]
        % Q           % Net heat flux at the water surface              [ W m-1 ]
        % Sn_0        % Shortwavera diation flux at the water surface   [ W m-1 ]
    %% OUTPUT %% ----------------------------------------------------------
        % T_new       % Water tempearture profile after the integation  [ C ]   


    %% ------------ LAYER CONFIFGURATION --------------------------------------
    % The depth of the integration points where water temperatures are calculated
    % are included in the input (z). Layer boundaries are set halfway between 
    % each two consecutive points. The layer number is calculated from the top
    % to the bottom.
        J  = length(z) ;    % Number of lake layers in the model
        dz = zeros(J,1);    % Thickness od j-th layer 
        zz = zeros(J,1);    % Depth of the boundary between j-th i j+1-th layer
        for j = 1:J
            if j==J
                zz(j)=z(j)+(z(j)-z(j-1))/2;
            else
                zz(j) = z(j)+(z(j+1)-z(j))/2;
            end

            if j==1
                dz(j)=zz(j);
            else
                dz(j)=zz(j)-zz(j-1);
            end
        end
    %% ------------ CONSTANTS -------------------------------------------------
        c_p  = 4188   ;     % Specific heat capacity of water           [ J kg-1 K-1]
        k_m  = 0.6    ;     % Molecular thermal conductivity of water	[ W m-1 K-1 ]
        g    = 9.8065 ;     % Gravitational acceleration                [ m s-2 ]
        k    = 0.4    ;     % Von Karman constant                       [ - ]
        Pr0  = 1      ;     % Neutral value of turbulent Prandtl number [ - ]
        beta = 0    ;     % Solar radiatio absorption at the surface  [ - ]
    %% ------------ VARIABLES -------------------------------------------------
        lambda_w = 1.1925*zz.^(-0.424);         % Light extinction coefficient              [ m-1 ] 
        Fi = (1-beta)*exp(-lambda_w.*zz)*Sn_0'; % Atenuated shortwave radiation at the 
                                                % boundary between j-th i j+1-th layer      [ W m-2 ]
        ro = ROw(T_old);                        % Water density in the integration points	[ kg m-3 ]        
        %% Turbulent diffusion 
        % Calculating the turbulent difusion thermal conductivity at each bondary between layers - k_t
            if TURB_DIF == 1        
                U2m = u_star/k*log(2/z0) ;                  % Wind speed at 2m hight above ground   [ m s-1 ]
                k_star = 6.6*sqrt(sind(LAT))*U2m^(-1.84);   % Parameter of the Ekman layer          [-]
                for j=1:J-1            
                    ro_j = (ro(j)+ro(j+1))/2;               % Water density at the layer boundaries	[ kg m-3 ]
                    d_ro = -(ro(j)-ro(j+1))/(z(j)-z(j+1));  % Water density gradient                [ kg m-4 ]
                    N2 = -g/ro_j*d_ro ;                     % Brunt-Vaisala frequency               [ Hz ]   
                    if N2>0 
                        Ri(j) = (-1 + sqrt((1+40*N2*k_m^2*zz(j)^2)/(u_star^2*exp(-2*k_star*zz(j)))))/20; 
                        k_t(j) = c_p * ro_j * ( k * u_star * zz(j)/Pr0 ) * exp(-k_star *zz(j)) / (1+37*Ri(j)^2);
                    else
                        k_t(j) = c_p * ro_j * ( k * u_star * zz(j)/Pr0 ) * exp(-k_star *zz(j)) ;
                    end
                end
            else
                k_t = zeros(J-1,1);
            end            
    %% --------- INTEGRATION STEP ---------------------------------------------
    %   Solving the equation  :    M * T(n+1) = A * T(n) - B 
        %% Matrix preparation
            A = diag(c_p*ro/dt.*dz); 
            M = zeros(J);
            B = zeros(J,1); 
            M(1,1) = c_p*ro(1)*dz(1)/dt + (k_m+k_t(1))/(z(2)-z(1)) ;  
            M(1,2) = -(k_m+k_t(1))/(z(2)-z(1));                
            B(1)   = Q - Fi(1) ;    
            for j = 2:J-1
                M(j,j-1) = -(k_m+k_t(j-1))/(z(j)-z(j-1)) ;                              
                M(j, j ) = c_p*ro(j)*dz(j)/dt + (k_m+k_t(j))/(z(j+1)-z(j))+(k_m+k_t(j-1))/(z(j)-z(j-1)) ; 
                M(j,j+1) = -(k_m+k_t(j))/(z(j+1)-z(j)) ;                               
                B(j  )   = Fi(j-1) - Fi(j)  ;    
            end
            M(J,J-1) = -(k_m+k_t(J-1))/(z(J)-z(J-1)) ;               
            M(J, J ) = c_p*ro(J)*dz(J)/dt + (k_m+k_t(J-1))/(z(J)-z(J-1)) ;  
            B(J)     = Fi(J-1)   ;    
        %% Integration 
             T_new = M \ (A * T_old) + M\B ;
        %% Convective mixing
            ConvMix = CONV_MIX;
            while ConvMix==1
                for j=2:J
                    ro_new = ROw(T_new);
                    d_ro = ro_new(j-1)-ro_new(j);
                    if d_ro > 0.0001
                            T_new(j-1) = (T_new(j-1)*dz(j-1) + T_new(j)*dz(j))/(dz(j-1)+dz(j));
                            T_new(j)   = T_new(j-1) ;
                            break
                    end
                    if j==J
                        ConvMix = 0 ;
                    end
                end
            end    
end
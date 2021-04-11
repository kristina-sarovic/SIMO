function [T_s] = SIMO_main(LAT, LON,z_ASL, Tmean, dTmean, T_initial, z, V,Ta,p,rh, UVB, P, f, DateTime, dt, i_s,i_f, TURB_DIF, CONV_MIX)

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
    % Main function ... 
    %
    %% INPUT %%
        % LAT       % Location latitude                 [ deg ]
        % LONG      % Location longitude                [ deg ]
        % z_ASL     % Loaction altitude                 [ m ]
        % Tmean     % Mean annual air temperature       [ C ]
        % dTmean    % Mean annual range between the... 
                    % ...daily air temp. max. and min.  [ C ]
        % T_initial % Initial water tempearture profile [ C ]   
        % z         % Depths of the integration points  [ m ]
        % V         % Wind speed                        [ m s-1 ]
        % Ta        % Air temerature                    [ C ] 
        % p         % Atmospheric pressure              [ Pa ]
        % rh        % Air relative humidity             [ - ]
        % UVB 	    % UVB radiaiton                     [ W m-2 ]
        % P         % Precipitation                     [ mm h-1 ]
        % f         % Cloud fraction                    [ - ]
        % DateTime  % Date and time         
        % dt        % Time step (dt=1h=3600s recomended)            [ s ]  
        % i_s       % Data row of the starting hour of the simulation 
        % i_f       % data row of the final hour of the simulation 
        % TURB_DIF  % Turbulent difusion flag (included/neglected) 	[ 1/0 ]
        % CONV_MIX  % Convective mixing  flag (included/neglected)  [ 1/0 ]
        
                
    %% OUTPUT %%
        % Tw_s      % Simulated water tempearture [ C ]  
        
    %% Data length check 
        N = length(V); % Simulation lenght [ h ]
        if mod(N,24*3600/dt)>0
            error( 'Input meteo data does not consist of full days.')
        end
        if ~(N==length(Ta) && N==length(p) && N==length(rh) && N==length(UVB) && N==length(P) && N==length(P)) 
            error( 'Input meteo data is not of consistent length.')
        end
        if ~(length(z)==length(T_initial)) 
            error( 'Integration point deapth vector and initial temperature vector are not of the same size.')
        end
                           
    %% Calculation 
        i  = i_s ;
        ii = 1 ;
        T_old = T_initial ; 
        T_s  = nan(i_f-i_s+1, length(T_initial)) ;  % Initialization of the result vector
        T_s(ii,:) = T_initial;
        while i<=i_f 
                temp = DateTime(i);
                dd_s = find(DateTime == datetime(year(temp), month(temp),day(temp), 0,0,0));   % first row of the day
                dd_k = find(DateTime == datetime(year(temp), month(temp),day(temp),23,0,0));   % last row of the day
            %% Calculating net surface heat fluxes    
            try
                [Q, Sn_0, ~, ~, ~, ~, u_star, z0] = SHF( LAT, LON, z_ASL, Tmean, dTmean,DateTime(i), dt, ...
                                                         V(i),  p(i), rh(i), T_old(1),  f(i),...
                                                         Ta(dd_s:dd_k), UVB(dd_s:dd_k), P(dd_s:dd_k) ) ;
            catch
                    a=1;
                end
            %% Calculating the water temperature in the next time step 
                 T_new = SIMO_compute_Tnew(T_old, z, dt, TURB_DIF, CONV_MIX, LAT, u_star, z0, Q, Sn_0)   ;             
            %% Saving the new water temperature and preparing for the next step
                i  = i+1 ;
                ii = ii+1;  
                T_old     = T_new ; 
                T_s(ii,:) = T_new' ;
        end
end
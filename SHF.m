function [Q, Sn_0, Ln, Hs, LE, Rp, u_star, z_0] = SHF( LAT, LON, z_ASL, Tmean, dTmean, DateTime, dt, V, p, rh, Ts, f, Ta_d, UVB_d, P_d)

%% Heat Exchange Model 
    %% Last update 26.03.2021 
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
    % This function calculates the surface heat flux components.
    %
    %% INPUT %% -----------------------------------------------------------
        % LAT       % Location latitude                 [ deg ]
        % LONG      % Location longitude                [ deg ]
        % z_ASL     % Loaction altitude                 [ m ]
        % Tmean     % Mean annual air temperature       [ C ]
        % dTmean    % Mean annual range between the... 
                    % ...daily air temp. max. and min.  [ C ]
        % DateTime  % Date and time
        % dt        % Time step (dt=1h=3600s recomended)[ s ]   
        % V         % Wind speed                        [ m s-1 ]
        % p         % Atmospheric pressure              [ Pa  ]   
        % rh        % Air relative humidity             [ - ]
        % Ts        % Surface water temperature         [ C ] 
        % f         % Cloud fraction                    [ - ]
        % Ta_d      % Air temperature (vector with values for the whole current day) [ C ]
        % UVB_d     % UVB radiation   (vector with values for the whole current day) [ W m-2 ] 
        % P_d       % Precipitation   (vector with values for the whole current day) [ mm s-1 ]
    %% OUTPUT %% ----------------------------------------------------------
        % Q         % Net surface heat flux         [ W m-2 ]
        % Sn_0      % Net shortwave radiation flux  [ W m-2 ]
        % Ln        % Net longwave radiation flux   [ W m-2 ]
        % Hs        % Sensible heat flux            [ W m-2 ]
        % LE        % Latent heat flux              [ W m-2 ]
        % Rp        % Precipitation heat flux       [ W m-2 ]
        % u_star    % Air shear velocity            [ m s-1 ]
        % z_0       % Roughness length              [ m ]
        
    %% --- Data check/edit ------------------------------------------------
        N = length(Ta_d);    
        if ~(N==24*3600/dt)
            error( 'Input meteo data does not represent a full day.')
        end
        if ~(N==length(UVB_d) && N==length(P_d))
            error( 'Daily input meteo data is not of consistent length.')
        end
        % If precipitations data includes NaNs for the periods without precipitation, convert them to zeros
            P_d(isnan(P_d))=0;
            
    %% --- CONSTANTS ------------------------------------------------------
        c_p     = 4188 ;    % Specific heat capacity of water [ J kg-1 K-1]
        E_w     = 0.96 ;    % Emissivity of water surface     [ - ]
        alpha_w = 0.06 ;    % Water surface albedo            [ - ]
        sigma =5.67*10^(-8);% Stefan-Boltzman constant        [ kg s-3 K-4 ]
  
    %% --- Calculation of variables ---------------------------------------
        n    = hour(DateTime)+1 ;       % Hour in the day for which the calulation is performed
        P    = P_d(n) ;                 % Precipitation in the hour of calculation      [ mm h-1 ]
        Psum = sum(P_d) ;               % Total daily precipitation                     [ mm day-1 ]
        Ta   = Ta_d(n) ;                % Air temperature in the hour of calculation	[ C ] 
        Tmax = nanmax(Ta_d) ;           % Daily air temperature maximum                 [ C ] 
        Tmin = nanmin(Ta_d(1:length(Ta_d)/2,:)); % Morning air temperature minimu       [ C ]
        ro_w = ROw(Ts) ;                         % Water density at the surface         [ kg m-3 ]
        es   = sat_pres(Ta);                     % Water vapour saturation pressure     [ Pa ]
        e    = rh*es ;                           % Water vapour partial pressure        [ Pa ]
        E_a  = 1.24*((e/100)/(Ta+273.15))^(1/7); % Emissivity of the atmosphere         [ - ]    
       
    %% --- CALCULATION OF THE HEAT FLUX COMPONENTS ------------------------
        %% Net shortwave radiation
            Sn_d   = solar_surf(LAT, LON, z_ASL, Tmean, dTmean, Psum, Tmin, Tmax, DateTime) ;  % Total daily solar irradiance at the surface   [ J day-1 m-2 ]
            UVB_wf = UVB_d / sum(UVB_d) ;               % Weight factors for calculating the hourly values of solar radiation
            Sn     = UVB_wf * Sn_d / 3600 ;             % Solar radiation reaching the water surface    [ W m-2 ]
            Sn_0   = Sn(n) * (1-alpha_w);           % Solar radiation absorbed at the water surface [ W m-2 ]
        %% Net longwave radiation            
            L_ac = E_a.*sigma.*(Ta+273.15).^4 ;     % Downward radiation from the atmosphere...  
                                                    %  ... under clear sky conditions           [ W m-2 ] 
            L_a = L_ac*(1-f)+f*sigma*(Ta+273.15)^4; % Downward radiation from the atmosphere    [ W m-2 ] 
            L_w  = E_w*sigma*(Ts+273.15)^4 ;      	% Upward radiation from the water surface   [ W m-2 ] 
            Ln   = E_w*L_a - L_w ;                  % Net longwave radiation                    [ W m-2 ]        
        %% Latent and Sensible heat flux
            temp    = sens_latent_heat(V,Ta,Ts,rh*100,p/100,10,2);          
            Hs      = -temp(:,14);      % Sensible heat flux 	[ W m-2 ] 
            LE      = -temp(:,12);      % Latent heat flux      [ W m-2 ]
            u_star  =  temp(:,5) ;      % Air shear velocity    [ m s-1 ]
            z_0     =  temp(:,15);      % Roughness length      [ m ]
        %% Heat flux brought by precipitation 
            % Assuming the precipitation temperature equals the air temperature
            Rp  = ro_w *P/1000/3600*Ta*c_p ;    % Precipitation heat flux 	[ W m-2 ]       
    %% NET TOTAL SURFACE HEAT FLUX [ W m-2 ]
        Q = Sn_0 + Ln + Hs + LE + Rp ; 

end
    
    
    
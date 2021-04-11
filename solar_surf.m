function Sn_d = solar_surf(LAT, LON, z_ASL, Tmean, dTmean, Psum, Tmin, Tmax, DateTime)

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
    % This function calculates the total daily solar irradiance at the Earth's surface
    %
    %% INPUT --------------------------------------------------------------
        % LAT       % Location latitude               [ deg ]
        % LON       % Location longitude              [ deg ]
        % z_ASL     % Loaction altitude               [ m ]
        % Tmean     % Mean annual air temperature     [ C ]
        % dTmean    % Mean annual range between the... 
                    % ...daily air temp. max. and min.[ C ]
        % DateTime  % Date and time
        % Psum      % Total daily precipitation       [ mm day-1 ]
        % Tmin      % Morning air temperature minimum [ C ] 
        % Tmax      % Daily air temperature maximum   [ C ]    
    %% OUTPUT -------------------------------------------------------------
        % Sn_d    	% Daily solar irradiance at earth's surface [ J day-1 m-2 ]
    
    %% CALCULATION --------------------------------------------------------   
        %% Daily solar irradiance at the top of the atmosphere and half day length
            [S_top, omega_s] = solar_top( DateTime, LAT, LON );
                % S_top     % Daily solar irradiance at the top of the atmosphere   [ J day-1 m-2 ]
                % omega_s   % Half day length                                       [ rad ]
        %% Daily total atmospheric transmitance to solar radiation  TAU [ - ]            
            % Transmiitance of dry clear air
                if abs(LAT)<80
                    tau_0 = 0.947 - (1.033 * 10^(-5))*abs(LAT)^2.22;
                else 
                    tau_0 = 0.774 ; 
                end
            % Transmiitance affected by aerosols (very variable, actually location dependant)
                tau_a = 1 ;                 
            % Transmiitance affected by water vapour
                tau_v = 0.9636 - 9.092 * 10^(-5) * (Tmean + 30)^1.8232 ;
                if Psum > 1
                    tau_v = tau_v - 0.13 ;
                end
            % Elevation correction 
                c_elev = (1 - 2.2569 * 10^(-5) * z_ASL )^5.2553  ;          
            % Total atmospheric transmitance    
                TAU = (tau_0 * tau_a * tau_v)^c_elev;      
        %% Day length correction  
            D = 1/(1-(omega_s-pi/4)^2/(2*omega_s^2)) ; 
        %% Site variation coefficient 
            beta = max( 1.041 , 23.753* dTmean / (Tmean+273.16) ) ; 
        %% Water vapour saturation pressure
            es_Tmin = sat_pres(Tmin) ; % Water vapour saturation pressure at the moment of T_min  [ Pa ]
            es_Tmax = sat_pres(Tmax) ; % Water vapour saturation pressure at the moment of T_max  [ Pa ]            
        %% Daily solar irradiance at earth's surface    [ J day-1 m-2 ]
            Sn_d = S_top * TAU * D * (1 - beta * (es_Tmin / es_Tmax )) ;
        %% Minimum value of Sn_d is set to 0.1*S_top
            if Sn_d < 0.1 * S_top
                Sn_d = 0.1 * S_top ;
            end
end
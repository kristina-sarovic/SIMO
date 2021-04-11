function [S_top, omega_s] = solar_top( DateTime, LAT, LON )

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
    % This function calculates the total daily solar irradiation at the top of 
    % the atmosphere and the Half day length.
    % Reference: Wald, L.: Basics in Solar radiation at Earth’s surface, hal-02164311, 2019.
    %
    %% INPUT %% -----------------------------------------------------------
        % LAT       % Location latitude               [ deg ]
        % LON       % Location longitude              [ deg ]
        % DateTime  % Date and time         
    %% OUTPUT %% ----------------------------------------------------------
        % Sn_d    	% Daily solar irradiance at earth's surface [ J day-1 m-2 ]
        % omega_s   % Half day length                           [ rad ]
        
    %% --- CONSTANTS ------------------------------------------------------
        S0m = 1362; % Solar constant  [ W m-2 ]
        
    %% CALCULATION %% -----------------------------------------------------
        % Check if DateTime is in datetime format
            if ~isdatetime(DateTime)
                error('Input DateTime is not in a datetime format')       
            end
        % Calculation of year, day of year, day angle and eccentricity factor
            y = year(DateTime);             % Year
            j = day(DateTime, 'dayofyear'); % Day of year
            jj  = j *2*pi/365.25 ;          % Day angle             [ rad ]
            ecc = 1+0.03344*cos(jj-0.049);  % Eccentircity factor
        % Convert locattion coordinates to radians
            lat = LAT*pi/180;   % Location latitude [ rad ]
            lon = LON*pi/180;   % Location latitude [ rad ]
        %% Solar declination at noon
            % For calculating the daily dolar radiation it is recommended
            % to take the solar declination value at noon. 
            n0 = 78.8944+0.2422*(y-1957)-floor((y-1957)/4); % Spring-equinox time expressed in days ...
                                                            % ...from the beginning of the year at LON=0
            t1 = -0.5 - lon/(2*pi) - n0 ;     % Spring equinox time correction for longitude
            omega_t = (2*pi/365.2422)*(j+t1); % Day angle counted from the spring equinox    [ rad ]  

            delta =   0.0064979  ...
                    + 0.4059059*sin(omega_t) + 0.0020054*sin(2*omega_t) - 0.002988*sin(3*omega_t) ...
                    - 0.0132296*cos(omega_t) + 0.0063809*cos(2*omega_t) + 0.0003508*cos(3*omega_t) ; % [ rad ]

        %% Half day length  [ rad ]        
            omega_s = acos(-tan(lat)*tan(delta));  
        %% % Daily solar irradiance at earth's surface [ J day-1 m-2 ]
            S_top   = S0m*ecc*(3600*24)/pi*(cos(lat)*cos(delta)*sin(omega_s)+omega_s*sin(lat)*sin(delta));  
        
end

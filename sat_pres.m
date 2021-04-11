function es = sat_pres(T)
    % This function calculates the water vapour saturation pressure as a 
    % function of the air temperature.
    %% --- INPUT ----------------------------------------------------------
        % T 	% Air temperature                   [ C ]
    %% --- OUTPUT ---------------------------------------------------------
        % es    % Water vapour saturation pressure  [ Pa ]
    %% --- CALCULATION ----------------------------------------------------        
        if T > 0
            es = 6.1094 * exp (17.269 * T/(T+237.7) ) * 100 ; 
        else
            es = 6.1094 * exp (21.875 * T/(T+265.3) ) * 100 ; 
        end
end
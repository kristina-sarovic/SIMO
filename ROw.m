function ro = ROw(T)
  % This function calculates the water density as a function of its temperature.
    %% --- INPUT ----------------------------------------------------------
        % T 	% Water temperature 	[ C ]
    %% --- OUTPUT ---------------------------------------------------------
        % ro    % Water density         [ kg m-3 ]
    %% --- CALCULATION ----------------------------------------------------        
    ro = (1-1.9549*10^(-5)*abs(T+273.15-277).^1.68)*10^3;
end
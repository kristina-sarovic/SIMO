function ro = ROw(T)
    % T     C       temperatura vode
    % ro    kg/m3   gustoca vode 
    ro = (1-1.9549*10^(-5)*abs(T+273.15-277).^1.68)*10^3; % kg/m3 -  gustoca vode 
end
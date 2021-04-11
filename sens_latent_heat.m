function mm=sens_latent_heat(U_Z,T,T0,Rh,p,z,z_t)

% USAGE
%
%    mm=compute_senslatentheat(Dates,U_Z,T,T0,Rh,p,z,z_t);
%
% Script to compute surface thermodynamic fluxes (sensible and latent) using the 
% bulk aerodynamic approach taking into account atmospheric stability 
% , including roughness lengths of momentum, vapor and temperature
%
% This set of scripts was constructed by Piet Verburg and Jason P. Antenucci
%
% Please reference the following paper if this script is of use to you.
% Verburg, P. and J.P. Antenucci. 2010. Persistent unstable atmospheric boundary layer
% enhances sensible and latent heat loss in a tropical great lake: Lake Tanganyika. 
% Journal of Geophysical Research - Atmospheres 115, D11109. doi:10.1029/2009JD012839
%
% Input
%   U_Z:    Wind Speed at the measurement height (m/s)
%   T:      Air Temperature at the measurement height (degrees C)
%   T0:     Surface Water Temperature (degrees C)
%   Rh:     Relative Humidity (%)
%   p :     Air Pressure (mb)
%   z:      Wind Sensor Height (m)
%   z_t:    Temperature and humidity Sensor Height (m)
%
% Output
%   mm=[Rh T T0 U_Z u_star zeta L C_H C_HN C_D C_DN E  Evap H  z_0 z_E z_T del_theta rho_a rho_w]; 
%
% Last update 2010-10-11
% SOURCE:  https://niwa.co.nz/our-services/software/heat-fluxes-from-lakes

%% MODIFICATIONS by K. Šaroviæ:
% - Converting the U_Z=0 to very small value (0.00001) to avoid numerical 
%   problems associated with U_Z=0
% - Modified saturated vapour pressure formula to keep consitancy with the
%   rest of the model
% - Introducing temperature and humidity sensor height different from the
%   wind sensor hight. 

%%
% Constants used in the calculations
const_SpecificHeatAir = 1005;           % Units : J kg-1 K-1
const_vonKarman = 0.4;                 % Units : none
const_Gravity = 9.81;                   % Units : m s-2
const_Charnock = 0.013;                 % Units : none

% Convert all input to columns
U_Z=U_Z(:);
U_Z(U_Z==0)=0.0001; 
T=T(:);
T0=T0(:);
Rh=Rh(:);
p=p(:);
z=z(:);

% Step 2c - Compute saturated vapour pressure at air temperature
e_s = sat_pres(T)/100; % Units : mbar   
% Step 2d - Compute vapour pressure
e_a = Rh.*e_s./100; % Units : mb
%%% End step 2

%% Step 3 - Compute other values used in flux calculations
% Step 3a - Compute specific humidity
q_z = 0.622.*e_a./p; % Units: kg kg-1
% Step 3b - Compute saturated vapour pressure at water temperature
e_sat = sat_pres(T0)/100 ; % Units : mb ##REF## 
% Step 3c - Compute humidity at saturation (Henderson-Sellers 1986 eqn 36)
q_s = 0.622.*e_sat./p; % Units: kg kg-1
% Step 3d - Compute latent heat of vaporisation
L_v = 2.501e6-2370*T0; % Units : J kg-1 ** EQUATION FROM PIET ##REF##
% Step 3e - Compute gas constant for moist air
R_a = 287*(1+0.608*q_z); % Units : J kg-1 K-1
% Step 3f - Compute air density
rho_a = 100*p./(R_a.*(T+273.16)); % Units : kg m-3
% Step 3g - Compute kinematic viscosity of air 
v = (1./rho_a).*(4.94e-8*T + 1.7184e-5); % Units : m2 s-1
% Step 3h - Compute virtual air temperature and virtual air-water temperature difference
T_v = (T+273.16).*(1+0.61*q_z); % Units - K
T_ov = (T0+273.16).*(1+0.61*q_s); % Units - K
del_theta = T_ov - T_v;
% Step 3h - Compute water density 
rho_w = ROw(T0);
%%% End step 3

%% Step 4 - Compute initial estimates of neutral transfer coefficients. 
%%%         This contains an iteration loop 
% Step 4a - Compute initial approximation to AIR shear velocity
u_star = U_Z.*sqrt(0.00104+0.0015./(1+exp((-U_Z+12.5)/1.56))); % Amorocho and DeVries, initialise ustar using U_Z
% Step 4b - Compute initial roughness length for momentum
z_0 = (const_Charnock*u_star.^2./const_Gravity) + (0.11*v./u_star); % Units: m
% Step 4c - Iterate initial u_star
% fprintf('Computing neutral coefficients\n');
z_0_prev=z_0*1.1; % To initiate the iteration
    for i1=1:length(U_Z)
        while abs((z_0(i1) - z_0_prev(i1)))/abs(z_0_prev(i1)) > 0.000001 % Converge when z_0 within 0.0001% of previous value
            u_star(i1)=const_vonKarman*U_Z(i1)/(log(z/z_0(i1)));  % Compute u_star
            dummy = z_0(i1); % Used to control while loop
            z_0(i1)=(const_Charnock*u_star(i1).^2./const_Gravity) + (0.11*v(i1)./u_star(i1)); % Compute new roughness length
            z_0_prev(i1) = dummy; % Used to control while loop
        end
    end
    % Step 4d - Compute initial neutral drag coefficient
C_DN = (u_star.^2)./(U_Z.^2); % Units - none
% Step 4e - Compute roughness Reynolds number 
Re_star = u_star.*z_0./v; % Units - none
% Step 4f - Compute initial roughness length for temperature
z_T = z_0.*exp(-2.67*(Re_star).^(1/4) + 2.57); % Units - m
z_T = real(z_T); % Get real components, and NaN can create imag component despite no data
% Step 4g - Compute initial roughness length for vapour 
z_E = z_0.*exp(-2.67*(Re_star).^(1/4) + 2.57); % Units - m
z_E = real(z_E); % Get real components, and NaN can create imag component despite no data
% Step 4h - Compute initial neutral sensible heat transfer coefficient 
C_HN = const_vonKarman*sqrt(C_DN)./(log(z_t./z_T));
% Step 4i - Compute initial neutral latent heat transfer coefficient
C_EN = const_vonKarman*sqrt(C_DN)./(log(z_t./z_E));
%%% End step 4

%%% Step 5 - Start iteration to compute corrections for atmospheric stability
% fprintf('Correcting neutral coefficients\n');
for i1=1:length(U_Z) % Need to iterate separately for each record
    % Step 5a - Compute initial sensible heat flux based on neutral coefficients
    H_initial(i1) = rho_a(i1)*const_SpecificHeatAir*C_HN(i1)*U_Z(i1)*(T0(i1)-T(i1)); % Units : W m-2
    % Step 5b - Compute initial latent heat flux based on neutral coefficients
    E_initial(i1) = rho_a(i1)*L_v(i1)*C_EN(i1)*U_Z(i1)*(q_s(i1)-q_z(i1)); % Units : W m-2
    % Step 5c - Compute initial Monin-Obukhov length
    L_initial(i1) = (-rho_a(i1)*u_star(i1)^3*T_v(i1))/(const_vonKarman*const_Gravity*(H_initial(i1)/const_SpecificHeatAir + 0.61*E_initial(i1)*(T(i1)+273.16)/L_v(i1))); % Units - m
    % Step 5d - Compute initial stability parameter
    zeta_initial(i1) = z./L_initial(i1);
    % Step 5e - Compute initial stability function
    psim(i1)=PSIM(zeta_initial(i1)); % Momentum stability function
    psit(i1)=PSITE(zeta_initial(i1)); % Sensible heat stability function
    psie(i1)=PSITE(zeta_initial(i1)); % Latent heat stability function
    % Step 5f - Compute corrected coefficients
    C_D(i1)=const_vonKarman*const_vonKarman/(log(z/z_0(i1))-psim(i1))^2;
    C_H(i1)=const_vonKarman*sqrt(C_D(i1))/(log(z_t/z_T(i1))-psit(i1));
    C_E(i1)=const_vonKarman*sqrt(C_D(i1))/(log(z_t/z_E(i1))-psie(i1));
    % Step 5g - Start iteration
    L_prev = L_initial(i1);
    L(i1) = L_prev*1.1; % Initialise while loop
    count(i1)=0;
     while abs((L(i1) - L_prev))./abs(L_prev) > 0.000001 % Converge when L within 0.0001% or previous L
            % Iteration counter
            count(i1)=count(i1)+1;
            if count(i1) > 20 ; break; end
            % Step 5i - Compute new z_O, roughness length for momentum
            z_0(i1)= (const_Charnock*u_star(i1).^2./const_Gravity) + (0.11*v(i1)./u_star(i1));
            % Step 5j - Compute new Re_star
            Re_star(i1) = u_star(i1).*z_0(i1)./v(i1);
            % Step 5k - Compute new z_T, roughness length for temperature
            z_T(i1) = z_0(i1)*exp(-2.67*(Re_star(i1)).^(1/4) + 2.57);
            % Step 5l - Compute new z_E, roughness length for vapour
            z_E(i1) = z_0(i1)*exp(-2.67*(Re_star(i1)).^(1/4) + 2.57);
            % Step 5p - Compute new stability parameter
            zeta(i1) = z./L(i1);
            %fprintf('zeta %g\n',zeta(i1));
            % Step 5q - Check and enforce bounds on zeta
            if zeta(i1) > 15; zeta(i1) = 15; elseif zeta(i1) < -15; zeta(i1) = -15; end
            % Step 5r - Compute new stability functions
            psim(i1)=PSIM(zeta(i1)); % Momentum stability function
            psit(i1)=PSITE(zeta(i1)); % Sensible heat stability function
            psie(i1)=PSITE(zeta(i1)); % Latent heat stability function
            % Step 5s - Compute corrected coefficients
            C_D(i1)=const_vonKarman*const_vonKarman/(log(z/z_0(i1))-psim(i1))^2;
            C_H(i1)=const_vonKarman*sqrt(C_D(i1))/(log(z_t/z_T(i1))-psit(i1));
            C_E(i1)=const_vonKarman*sqrt(C_D(i1))/(log(z_t/z_E(i1))-psie(i1));
            % Step 5m - Compute new H (now using corrected coefficients)
            H(i1) = rho_a(i1)*const_SpecificHeatAir*C_H(i1)*U_Z(i1)*(T0(i1)-T(i1));
            % Step 5n - Compute new E (now using corrected coefficients)
            E(i1) = rho_a(i1)*L_v(i1)*C_E(i1)*U_Z(i1)*(q_s(i1)-q_z(i1));
            % Step 5h - Compute new u_star
            u_star(i1)=sqrt(C_D(i1)*U_Z(i1)^2);
            % Step 5o - Compute new Monin-Obukhov length
            dummy = L(i1); % Used to control while loop
            L(i1) = (-rho_a(i1)*u_star(i1)^3*T_v(i1))/(const_vonKarman*const_Gravity*(H(i1)/const_SpecificHeatAir + 0.61*E(i1)*(T(i1)+273.16)/L_v(i1)));
            L_prev = dummy; % Used to control while loop
     end
end
%%% End step 5

% Take real values to remove any complex values that arise from missing data or NaN.
C_D=real(C_D);
C_E=real(C_E);
C_H=real(C_H);
z_0=real(z_0);
z_E=real(z_E);
z_T=real(z_T);

% Compute evaporation [mm/day]
Evap = 86400*1000*E(:)./(rho_w.*L_v);

% Clean up output signal based on any missing data points
index=find(isnan(Rh)); E(index)=NaN; H(index)=NaN; Evap(index)=NaN; zeta(index)=NaN; ustar(index)=NaN; C_H(index)=NaN; C_D(index)=NaN; C_HN(index)=NaN; C_DN(index)=NaN; L(index)=NaN; z_0(index)=NaN; z_E(index)=NaN; z_T(index)=NaN; del_theta(index)=NaN;
index=find(isnan(T)); E(index)=NaN; H(index)=NaN; Evap(index)=NaN; zeta(index)=NaN; ustar(index)=NaN; C_H(index)=NaN; C_D(index)=NaN; C_HN(index)=NaN; C_DN(index)=NaN; L(index)=NaN; z_0(index)=NaN; z_E(index)=NaN; z_T(index)=NaN; del_theta(index)=NaN;
index=find(isnan(T0)); E(index)=NaN; H(index)=NaN; Evap(index)=NaN; zeta(index)=NaN; ustar(index)=NaN; C_H(index)=NaN; C_D(index)=NaN; C_HN(index)=NaN; C_DN(index)=NaN; L(index)=NaN; z_0(index)=NaN; z_E(index)=NaN; z_T(index)=NaN; del_theta(index)=NaN;
index=find(isnan(U_Z)); E(index)=NaN; H(index)=NaN; Evap(index)=NaN; zeta(index)=NaN; ustar(index)=NaN; C_H(index)=NaN; C_D(index)=NaN; C_HN(index)=NaN; C_DN(index)=NaN; L(index)=NaN; z_0(index)=NaN; z_E(index)=NaN; z_T(index)=NaN; del_theta(index)=NaN;
index=find(isnan(p)); E(index)=NaN; H(index)=NaN; Evap(index)=NaN; zeta(index)=NaN; ustar(index)=NaN; C_H(index)=NaN; C_D(index)=NaN; C_HN(index)=NaN; C_DN(index)=NaN; L(index)=NaN; z_0(index)=NaN; z_E(index)=NaN; z_T(index)=NaN; del_theta(index)=NaN;

mm=[Rh T T0 U_Z u_star zeta(:) L(:) C_H(:) C_HN C_D(:) C_DN E(:) Evap(:) H(:) z_0 z_E z_T del_theta rho_a rho_w]; 
%   1  2 3  4   5      6       7    8      9    10     11   12   13      14   15  16  17  18        19    20     


%% Function list, used by script above
%-------------------------------------------------------
function psim=PSIM(zeta)
% Function to compute stability functions for momentum
   if zeta < 0.0
      X = (1 - 16*zeta)^0.25; 
      psim = 2*log((1 + X)/2) + log((1 + X*X)/2)-2*atan(X) + pi/2;      
   
   elseif zeta > 0.0   % Stable case
      if zeta > 0.5      
        if zeta > 10.0
           psim = log(zeta) - 0.76*zeta - 12.093;
        else
          psim = 0.5/(zeta*zeta) - 4.25/zeta - 7.0*log(zeta) - 0.852; 
        end
      else 
        psim = -5*zeta ;
      end
   else 
     psim = 0.0;
   end
return  
   
%-------------------------------------------------------
function psite=PSITE(zeta)
% Function to compute stability functions for sensible and latent heat
   if zeta < 0.0
      X = (1 - 16*zeta)^0.25; 
      psite = 2*log((1 + X*X)/2);
   
   elseif zeta > 0.0   % Stable case
      if zeta > 0.5      
        if zeta > 10.0
           psite = log(zeta) - 0.76*zeta - 12.093;
        else
          psite = 0.5/(zeta*zeta) - 4.25/zeta - 7.0*log(zeta) - 0.852; 
        end
      else 
        psite = -5*zeta ;
      end
   else 
     psite = 0.0;
   end
return

%-------------------------------------------------------

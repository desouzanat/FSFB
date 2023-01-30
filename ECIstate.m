function [rECI,vECI] = ECIstate(H, e, nu, Q_peri2ECI)

mu = 398600;   % 


rperi = (H^2/mu) * (1/(1+ (e * cosd(nu)))) .* [cosd(nu); sind(nu); 0]; % position in perifocal frame 
vperi = (mu/H) .* [-sind(nu); (e + cosd(nu)); 0]; % velocity in perifocal frame

rECI = Q_peri2ECI' *  rperi;   % position in ECI (from perifocal)
vECI = Q_peri2ECI' * vperi;    % velocity in ECI (from perifocal)


end
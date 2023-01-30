function [Q_peri2ECI] = EulerRotation313(RAAN, inc, AOP)

% this functions does a 3-1-3/z-x-z/AOP-RAAN-AOP rotation to go from the
% perifocal to ECI frames

Q_peri2ECI = zeros(3,3);

Q_peri2ECI(1,1) = (cosd(AOP) * cosd(RAAN)) - (sind(RAAN) * cosd(inc) * sind(AOP));
Q_peri2ECI(1,2) = (cosd(AOP) * sind(RAAN)) + (cosd(RAAN) * cosd(inc) * sind(AOP));
Q_peri2ECI(1,3) = sind(AOP) * sind(inc);
Q_peri2ECI(2,1) = (-sind(AOP) * cosd(RAAN)) - (sind(RAAN) * cosd(AOP) * cosd(inc));
Q_peri2ECI(2,2) = (-sind(AOP) * sind(RAAN)) + (cosd(RAAN) * cosd(AOP) * cosd(inc));
Q_peri2ECI(2,3) = cosd(AOP) * sind(inc);
Q_peri2ECI(3,1) = sind(inc) * sind(RAAN);
Q_peri2ECI(3,2) = -sind(inc) * cosd(RAAN);
Q_peri2ECI(3,3) = cosd(inc);


end
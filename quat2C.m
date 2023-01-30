function [C] = quat2C(q)

% this function calculates the corresponding direction cosine matrix given
% a quaternion q

eta = q(4);

epsX = [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0];
eps = [q(1) q(2) q(3)]';

C = (((2 * eta^2) - 1) * eye(3)) + (2 * (eps) * eps') - (2 * eta * epsX);
end

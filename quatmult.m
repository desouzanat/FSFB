function [q3] = quatmult(q1,q2)

% this function computes the quaternion multiplication product between two
% quaternions

eps1 = q1(1:3);
eta1 = q1(4);

eps2 = q2(1:3);
eta2 = q2(4);

q3 = zeros(4,1);

eps1X = [0 -q1(3) q1(2); q1(3) 0 -q1(1); -q1(2) q1(1) 0];

q3(1:3) = (eta1 * eps2) + (eta2 * eps1) + (eps1X * eps2);
q3(4) = (eta1 * eta2) - (eps1' * eps2);

end

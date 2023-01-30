function qdot = quatdiff(q, w)

% quaternion differentiation equation
eta = q(4);
eps = q(1:3);
epsX = [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0];

etadot = -eps' * w / 2;
epsdot = (eta * eye(3) + epsX) * w / 2;

qdot = [epsdot; etadot];

end

function qstar = quatconj(q)

% quaternion conjugate equation
eta = q(4);
eps = q(1:3);
etastar = eta;
epsstar = -eps;

qstar = [epsstar; etastar];

end
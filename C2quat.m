function [quat] = C2quat(C)

% this function calculates the corresponding attittude quaternion given a
% direction cosine matrix C

eta = (1/2) * sqrt(1 + trace(C));
eps1 = (1/4) * ((C(2,3)- C(3,2))/eta);
eps2 = (1/4) * ((C(3,1)- C(1,3))/eta);
eps3 = (1/4) * ((C(1,2)- C(2,1))/eta);

quat = [eps1 eps2 eps3 eta]';

end
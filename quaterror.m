function qe = quaterror(q_b_a, qcommand_b_a, q_c_b)

% this function calculates the error quaternion given three initial
% quaternions

q_c_a = quatmult(quatconj(q_c_b), q_b_a);
qcommand_c_b= quatmult(quatconj(q_c_b), qcommand_b_a); % target quaternion

qe = quatmult(quatconj(qcommand_c_b), q_c_a);

end
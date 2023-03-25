function w = log_map(R)
phi = acos((trace(R) - 1)/2);

log_R = phi / 2*sin(phi) * (R - transpose(R));
w = [-log_R(2,3), log_R(1,3), -log_R(1,2)]';
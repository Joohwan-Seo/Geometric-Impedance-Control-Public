function err = err_fun(g_,gd_)
%     psi = eye(4) - inv(reshape(gd_,[4,4])) * reshape(g_,[4,4]);
%     err = trace(psi'*psi);
%     err = norm(psi,'fro');

    g = reshape(g_,[4,4]); gd = reshape(gd_,[4,4]);

    err = trace(eye(3) - gd(1:3,1:3)'*g(1:3,1:3)) + 1/2 * (g(1:3,4) - gd(1:3,4))'*(g(1:3,4) - gd(1:3,4));
end
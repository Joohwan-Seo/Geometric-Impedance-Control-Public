function M_mat = M_fun(q1,q2,q3,q4,q5,q6)
%M_FUN
%    M_MAT = M_FUN(Q1,Q2,Q3,Q4,Q5,Q6)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    01-Oct-2022 21:26:48

t2 = cos(q3);
t3 = cos(q4);
t4 = cos(q5);
t5 = sin(q2);
t6 = sin(q3);
t7 = sin(q4);
t8 = sin(q5);
t9 = q2+q5;
t10 = q3+q4;
t11 = q2.*2.0;
t12 = q3.*2.0;
t13 = q4.*2.0;
t14 = q5.*2.0;
t20 = -q5;
t40 = atan(3.818188629788691e-2);
t42 = atan(7.656026541520085e+2);
t44 = 9.934139466657452e+7;
t47 = 1.077033123737933e+9;
t15 = t4.^2;
t16 = q2+t10;
t17 = q3+t9;
t18 = sin(t9);
t19 = sin(t10);
t21 = -t14;
t23 = t9+t10;
t27 = q2+t20;
t37 = t4.*3.38e-2;
t41 = t6.*4.0271725e-2;
t45 = -t42;
t46 = t5.*3.53883322e-1;
t48 = t2.*1.054733773125;
t50 = t2.*t7.*6.92336305e-2;
t51 = t3.*t6.*6.92336305e-2;
t55 = t7.*1.2779713677e-1;
t56 = t7.*6.3898568385e-2;
t60 = t7.*7.747599173958887e+27;
t64 = t6.*t7.*t8.*1.45620788325e-2;
t65 = t4.*3.462224718e-3;
t66 = t2.*t3.*t8.*1.45620788325e-2;
t71 = t3.*t8.*1.3439942169525e-2;
t72 = t3.*t8.*2.687988433905e-2;
t22 = cos(t16);
t24 = sin(t16);
t25 = sin(t17);
t26 = cos(t23);
t28 = q3+t27;
t29 = sin(t27);
t30 = q5+t23;
t31 = t16+t20;
t32 = t16+t21;
t43 = -t41;
t49 = q2+q3+t45;
t52 = -t50;
t53 = -t51;
t57 = -t55;
t58 = -t56;
t61 = t19.*8.394466918884709e+27;
t62 = t18.*7.28103941625e-3;
t67 = -t64;
t68 = -t65;
t70 = t60-1.995836665581788e+27;
t86 = t15.*1.52851393612391e-2;
t33 = sin(t28);
t34 = cos(t30);
t35 = cos(t31);
t36 = cos(t32);
t38 = t8.*t24.*3.38e-2;
t54 = cos(t49);
t59 = t22.*2.1793103212e-2;
t63 = -t62;
t69 = t29.*7.28103941625e-3;
t73 = t25.*6.7199710847625e-3;
t75 = t26.*4.0267812573e-3;
t76 = t26.*5.645565393e-4;
t84 = t61+t70-1.995836665581788e+27;
t85 = (t4.*t70)./5.764607523034235e+29;
t87 = t22.*2.379148606387609e-1;
t93 = t58+t71+t86+3.639621178351609e-1;
t97 = t43+t48+t52+t53+t57+t66+t67+t72+t86+1.531446074036411;
t39 = -t38;
t74 = -t73;
t77 = t47.*t54.*1.3e-10;
t78 = t33.*6.7199710847625e-3;
t79 = -t75;
t80 = -t76;
t81 = t35.*4.0267812573e-3;
t82 = t35.*5.645565393e-4;
t88 = -t87;
t89 = t34.*3.821284840309775e-3;
t90 = t36.*3.821284840309775e-3;
t92 = (t4.*t84)./5.764607523034235e+29;
t95 = t52+t53+t66+t67+t93;
t83 = -t82;
t91 = -t90;
t94 = t59+t80+t81+t89+t91;
t96 = t63+t69+t74+t78+t79+t83+t88;
t98 = t73+t77+t78+t94;
t99 = t46+t62+t69+t98;
et1 = t4.*9.1826755932e-3-t19.*6.92336305e-2+t58+cos(t10.*2.0+t11).*4.285508656149022e-2;
et2 = sin(q3+t10+t11).*(-6.3898568385e-2)-sin(t10+t11+t20).*7.28103941625e-3+cos(t23.*2.0).*1.910642420154888e-3;
et3 = cos(t10.*2.0+t11+t20).*(-1.731112359e-3)+cos(t10.*2.0+t11+t21).*1.910642420154888e-3;
et4 = sin(q3+t10+t11+t20).*(-6.7199710847625e-3)+cos(t16+t23).*1.731112359e-3+sin(q4+q5).*6.7199710847625e-3;
et5 = sin(q5+t10).*7.28103941625e-3-sin(q4+t20).*6.7199710847625e-3-sin(t10+t11).*6.92336305e-2;
et6 = sin(t9+t16).*7.28103941625e-3-sin(t10+t20).*7.28103941625e-3+sin(t16+t17).*6.7199710847625e-3+cos(t11).*1.195796359375;
et7 = cos(t14).*(-3.821284840309775e-3);
et8 = (3.267693522137215e+25.*cos(t11+t12+atan(6.416451315410032e-2)))./5.62949953421312e+25+t44.*cos(q3+t40).*1.0625e-8+t44.*cos(q3+t11+t40).*1.0625e-8;
et9 = 2.264011855385696;
et10 = t2.*2.10946754625-t6.*8.054345e-2+t57+t72+t86-t2.*t7.*1.38467261e-1-t3.*t6.*1.38467261e-1+t2.*t3.*t8.*2.9124157665e-2;
et11 = t6.*t7.*t8.*(-2.9124157665e-2)+3.938138792786411;
mt1 = [et1+et2+et3+et4+et5+et6+et7+et8+et9,t99,t98,t94,t96,t39,t99,et10+et11,t97,t95,t92,t37,t98,t97,t57+t72+t86+1.531446074036411,t93,t85,t37,t94,t95,t93];
mt2 = [t8.^2.*(-1.52851393612391e-2)+3.792472571964e-1,t68,t37,t96,t92,t85,t68];
mt3 = [2.379148606387609e-1,0.0,t39,t37,t37,t37,0.0,3.38e-2];
M_mat = reshape([mt1,mt2,mt3],6,6);

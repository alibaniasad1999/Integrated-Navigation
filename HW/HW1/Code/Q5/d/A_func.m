function out1 = A_func(in1,in2)
%A_func
%    OUT1 = A_func(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    20-Dec-2023 16:11:47

L = in1(1,:);
a_x = in2(1,:);
a_y = in2(2,:);
a_z = in2(3,:);
h = in1(3,:);
phi = in1(7,:);
psi = in1(9,:);
theta = in1(8,:);
v_d = in1(6,:);
v_e = in1(5,:);
v_n = in1(4,:);
w_y = in2(5,:);
w_z = in2(6,:);
t2 = cos(L);
t3 = sin(L);
t4 = tan(L);
t5 = cos(phi);
t6 = cos(psi);
t7 = cos(theta);
t8 = sin(phi);
t9 = sin(psi);
t10 = sin(theta);
t11 = tan(theta);
t12 = L.*2.0;
t13 = h.^2;
t14 = h.^3;
t16 = v_e.^2;
t17 = v_n.^2;
t39 = sqrt(2.0);
t40 = L./2.0;
t41 = psi./2.0;
t42 = theta./2.0;
t66 = h.*1.073741824e+9;
t72 = 2.317047500592079e+4;
t78 = 1.472899174203877e+11;
t15 = t13.^2;
t18 = t2.^2;
t19 = t2.^3;
t21 = t2.^5;
t23 = sin(t12);
t24 = t3.^2;
t25 = t4.^2;
t26 = t5.^2;
t27 = t8.^2;
t28 = a_x.*t10;
t29 = t5.*w_y;
t30 = t8.*w_z;
t31 = t5.*t6;
t32 = 1.0./t2;
t34 = t5.*t9;
t35 = t6.*t8;
t36 = t8.*t9;
t37 = 1.0./t7;
t43 = sin(t40);
t44 = sin(t41);
t45 = sin(t42);
t46 = a_z.*t5.*t7;
t47 = a_y.*t7.*t8;
t80 = t2.*1.458423302905242e-4;
t81 = t3.*1.458423302905242e-4;
t83 = t13.*3.278727125438313e+35;
t84 = t2.*7.29211651452621e-5;
t85 = t3.*7.29211651452621e-5;
t20 = t18.^2;
t22 = t18.^3;
t38 = t37.^2;
t48 = t25+1.0;
t49 = -t28;
t50 = -t30;
t51 = t10.*t36;
t52 = t18-1.0;
t53 = t10.*t31;
t54 = t10.*t34;
t55 = t10.*t35;
t56 = t43.^2;
t57 = t44.^2;
t58 = t45.^2;
t59 = t26+t27;
t76 = t18.*7.718094650639642e+15;
t77 = t18.*3.859047325319821e+15;
t87 = t24.*6.6943799901378e-3;
t88 = t13.*t18.*4.41939415613864e+33;
t112 = t24.*5.185915816319999e-2;
t113 = t23.*5.77038762e-5;
t33 = t20.^2;
t60 = -t54;
t61 = -t55;
t62 = t57.*2.0;
t63 = t58.*2.0;
t64 = t56.*2.0;
t65 = t29+t50;
t70 = t31+t51;
t71 = t36+t53;
t75 = t46+t47+t49;
t79 = t13.*t20.*1.489224625905806e+31;
t82 = t77+5.726017049781037e+17;
t86 = t76+1.145203409956207e+18;
t89 = -t87;
t95 = t87-1.0;
t114 = -t113;
t67 = t64-1.0;
t68 = t62-1.0;
t69 = t63-1.0;
t73 = t34+t61;
t74 = t35+t60;
t90 = sqrt(t82);
t94 = t89+1.0;
t96 = t86.^(3.0./2.0);
t99 = 1.0./t95;
t126 = (h.*t72.*t78.*t95)./2.169431977370464e+22;
t130 = t112+t114+9.780317999999999;
t91 = t90.^3;
t92 = t90.^5;
t93 = t90.^7;
t100 = t99.^3;
t101 = h.*t96;
t102 = 1.0./sqrt(t94);
t122 = h.*t39.*t90.*1.960736106370943e+33;
t127 = t126-1.0;
t97 = 1.0./t91;
t103 = t102.^3;
t104 = t102.^5;
t105 = t102.*6.378178e+6;
t111 = h.*t39.*t91.*3.424258239758336e+15;
t121 = t101+7.842944425483774e+33;
t128 = 1.0./t127.^3;
t98 = 1.0./t97.^3;
t106 = h+t105;
t110 = t103.*6.802669907834067e+15;
t115 = t103.*6.335480052823263e+6;
t123 = 1.0./t121;
t135 = t79+t83+t88+t111+t122+1.342813353646476e+49;
t107 = 1.0./t106;
t116 = h+t115;
t117 = t66+t110;
t136 = 1.0./t135;
t108 = t107.^2;
t109 = t107.*v_e;
t118 = 1.0./t116;
t120 = 1.0./t117.^2;
t119 = t118.^2;
t124 = t48.*t118.*v_e;
t125 = t84+t109;
t129 = t2.*t3.*t103.*t108.*v_e.*4.269794717673713e+4;
t131 = t2.*t3.*t4.*t104.*t119.*v_e.*1.272363326806117e+5;
t132 = t81+t129;
t134 = t85+t129;
t133 = -t131;
t137 = t80+t124+t133;
t138 = t84+t124+t133;
et1 = t132.*v_e-1.0./t127.^2.*(cos(t12).*1.154077524e-4-t2.*t3.*1.037183163264e-1);
et2 = t2.*t3.*t17.*t104.*t119.*1.272363326806117e+5-h.*t2.*t3.*t72.*t78.*t128.*t130.*1.234310189942339e-24;
et3 = t90.*1.803147702731296e+98;
et4 = h.*t39.*6.030418132156747e+100+t14.*t39.*1.472438105708829e+87;
et5 = t13.*t91.*7.688972157653375e+66+t13.*t92.*5.371253414585904e+49+t13.*t93.*2.345108898510572e+31+t15.*t98;
et6 = h.*t18.*t39.*6.096297503432047e+98;
et7 = h.*t20.*t39.*1.369532572143551e+96+t14.*t18.*t39.*3.473222135818082e+85;
et8 = t14.*t20.*t39.*3.009570334654626e+83+t14.*t22.*t39.*1.126832748113663e+81;
et9 = t14.*t33.*t39.*1.518857126999683e+78;
et10 = t6.*t90.*9.702074742074123e+113;
et11 = h.*t6.*t39.*3.244746248769314e+116;
et12 = t6.*t14.*t39.*7.922648007717839e+102+t6.*t13.*t91.*4.137153182198113e+82;
et13 = t6.*t13.*t92.*2.890073952788051e+65+t6.*t13.*t93.*1.261816864874061e+47+t6.*t15.*t98.*5.380632283965437e+15;
et14 = h.*t6.*t18.*t39.*3.280193515962437e+114;
et15 = h.*t6.*t20.*t39.*7.368951171617812e+111;
et16 = t6.*t14.*t18.*t39.*1.868813115336616e+101;
et17 = t6.*t14.*t20.*t39.*1.619339130350735e+99;
et18 = t6.*t14.*t22.*t39.*6.063072663129866e+96;
et19 = t6.*t14.*t33.*t39.*8.172411692265486e+93;
et20 = t2.*t6.*t39.*v_e.*7.497122656200364e+117;
et21 = t2.*t9.*t39.*v_n.*-2.249136796860109e+118;
et22 = t9.*t19.*t39.*v_n.*-1.515805011536474e+116;
et23 = t6.*t13.*1.0./t32.^7.*t39.*v_e.*5.603593741258584e+97+h.*t2.*t6.*t91.*v_e.*7.647253122784097e+84;
et24 = h.*t2.*t9.*t91.*v_n.*-2.294175936835229e+85;
et25 = t2.*t6.*t13.*t39.*v_e.*1.830561138587972e+104;
et26 = t6.*t13.*t19.*t39.*v_e.*3.701118249013487e+102;
et27 = t6.*t13.*t21.*t39.*v_e.*2.494367438199306e+100;
et28 = t2.*t9.*t13.*t39.*v_n.*-5.491683415763917e+104;
et29 = t9.*t13.*t19.*t39.*v_n.*-7.402236498026973e+102;
et30 = t9.*t13.*t21.*t39.*v_n.*-2.494367438199306e+100;
et31 = t8.*w_y.*1.342813353646476e+49+t5.*w_z.*1.342813353646476e+49-t2.*t6.*t10.*9.791951432051792e+44+t8.*t79.*w_y+t8.*t83.*w_y+t8.*t88.*w_y+t8.*t111.*w_y+t8.*t122.*w_y+t5.*t79.*w_z+t5.*t83.*w_z+t5.*t88.*w_z+t5.*t111.*w_z+t5.*t122.*w_z;
et32 = h.*t6.*t10.*v_e.*-3.278727125438313e+35+h.*t9.*t10.*v_n.*3.278727125438313e+35-t2.*t6.*t10.*t13.*2.390886021803377e+31;
et33 = t6.*t10.*t13.*t19.*(-3.222673711017921e+29)-t6.*t10.*t13.*t21.*1.085959948840685e+27-h.*t6.*t10.*t18.*v_e.*4.41939415613864e+33;
et34 = h.*t6.*t10.*t20.*v_e.*-1.489224625905806e+31+h.*t9.*t10.*t18.*v_n.*4.41939415613864e+33+h.*t9.*t10.*t20.*v_n.*1.489224625905806e+31-t6.*t10.*t39.*t90.*v_e.*1.960736106370943e+33+t9.*t10.*t39.*t91.*v_n.*3.424258239758336e+15;
et35 = h.*t2.*t6.*t10.*t39.*t90.*(-1.429791614189538e+29)-h.*t2.*t6.*t10.*t39.*t91.*2.497009006014421e+11;
et36 = t2.*t9.*9.791951432051792e+44+h.*t9.*v_e.*3.278727125438313e+35+h.*t6.*v_n.*3.278727125438313e+35;
et37 = t2.*t9.*t13.*2.390886021803377e+31+t9.*t13.*t19.*3.222673711017921e+29;
et38 = t9.*t13.*t21.*1.085959948840685e+27+h.*t9.*t18.*v_e.*4.41939415613864e+33+h.*t9.*t20.*v_e.*1.489224625905806e+31+h.*t6.*t18.*v_n.*4.41939415613864e+33;
et39 = h.*t6.*t20.*v_n.*1.489224625905806e+31+t9.*t39.*t90.*v_e.*1.960736106370943e+33+t6.*t39.*t91.*v_n.*3.424258239758336e+15+h.*t2.*t9.*t39.*t90.*1.429791614189538e+29+h.*t2.*t9.*t39.*t91.*2.497009006014421e+11;
et40 = t18.*3.79582894316288e+51+t20.*2.55819767966713e+49+t22.*5.746988309402213e+46;
et41 = 1.877404742183935e+53;
et42 = (t2.*t3.*t6.*t39.*v_n)./t92;
et43 = 1.0./(t13.*6.189700196426901e+26+4.759243247454807e+93./(et40+et41)+h.*t39.*t97.*2.427273732549109e+60);
et44 = -1.405044630811905e+76;
et45 = t3.*t9.*(-7.29211651452621e-5)+et42.*et43.*et44;
et46 = (t2.*t3.*t9.*t39.*t97.*v_e.*-3.808784161954939e+48)./(t13.*2.882303761517117e+17+1.0./t90.^2.*6.759316199344284e+48+(h.*t39.*1.973950480972287e+33)./t90);
et47 = t4.*t7.*t13.*v_e.*1.152921504606847e+18-t6.*t10.*t13.*v_e.*1.152921504606847e+18-t4.*t7.*t99.*v_e.*4.690217797021143e+31+t6.*t10.*t100.*v_e.*4.627631787495115e+31+t9.*t10.*t13.*v_n.*1.152921504606847e+18-t9.*t10.*t99.*v_n.*4.690217797021143e+31;
et48 = h.*t4.*t7.*t102.*v_e.*1.470707715282058e+25-h.*t6.*t10.*t103.*v_e.*1.460862238981533e+25+h.*t9.*t10.*t102.*v_n.*1.470707715282058e+25;
et49 = t2.*t6.*-9.791951432051792e+44-h.*t6.*v_e.*3.278727125438313e+35+h.*t9.*v_n.*3.278727125438313e+35;
et50 = t2.*t6.*t13.*(-2.390886021803377e+31)-t6.*t13.*t19.*3.222673711017921e+29;
et51 = t6.*t13.*t21.*(-1.085959948840685e+27)+t8.*t10.*w_y.*1.342813353646476e+49+t5.*t10.*w_z.*1.342813353646476e+49;
et52 = h.*t6.*t18.*v_e.*-4.41939415613864e+33-h.*t6.*t20.*v_e.*1.489224625905806e+31+h.*t9.*t18.*v_n.*4.41939415613864e+33+h.*t9.*t20.*v_n.*1.489224625905806e+31-t6.*t39.*t90.*v_e.*1.960736106370943e+33;
et53 = t9.*t39.*t91.*v_n.*3.424258239758336e+15+t8.*t10.*t79.*w_y+t8.*t10.*t83.*w_y+t8.*t10.*t88.*w_y+t8.*t10.*t111.*w_y+t8.*t10.*t122.*w_y+t5.*t10.*t79.*w_z+t5.*t10.*t83.*w_z+t5.*t10.*t88.*w_z+t5.*t10.*t111.*w_z+t5.*t10.*t122.*w_z-h.*t2.*t6.*t39.*t90.*1.429791614189538e+29-h.*t2.*t6.*t39.*t91.*2.497009006014421e+11;
mt1 = [t23.*t39.*t90.*v_n.*1.0./(t39.*1.960736106370943e+33+h.*t91).^2.*-1.134986014042318e+49,-(t3.*t109)./(t24-1.0)-t3.*t103.*t108.*v_e.*4.269794717673713e+4,0.0];
mt2 = [-t137.*v_e-t2.*t3.*t104.*t119.*v_d.*v_n.*1.272363326806117e+5,-t132.*v_d+t137.*v_n,et1+et2,(t3.*t37.*(et10+et11+et12+et13+et14+et15+et16+et17+et18+et19+et20+et21+et22+et23+et24+et25+et26+et27+et28+et29+et30))./(et3.*7.378697629483821e+19+et4.*7.378697629483821e+19+et5.*7.378697629483821e+19+et6.*7.378697629483821e+19+et7.*7.378697629483821e+19+et8.*7.378697629483821e+19+et9.*7.378697629483821e+19),et45+et46];
mt3 = [t5.*t37.*(t71.*t134+t5.*t7.*t138+t2.*t3.*t74.*t104.*t119.*v_n.*1.272363326806117e+5)-t8.*t37.*(t73.*t134-t7.*t8.*t138+t2.*t3.*t70.*t104.*t119.*v_n.*1.272363326806117e+5),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t119.*v_n,-t32.*t108.*v_e,0.0,t4.*t16.*t119-t119.*v_d.*v_n];
mt4 = [-t108.*v_d.*v_e-t4.*t119.*v_e.*v_n,t16.*t108+t17.*t119-t72.*t78.*t95.*t128.*t130.*9.219003042557575e-23,t5.*t11.*(t71.*t108.*v_e+t74.*t119.*v_n-t4.*t5.*t7.*t119.*v_e)-t8.*t11.*(t73.*t108.*v_e+t70.*t119.*v_n+t4.*t7.*t8.*t119.*v_e)+t6.*t7.*t108.*v_e+t4.*t10.*t119.*v_e-t7.*t9.*t119.*v_n];
mt5 = [-t59.*t108.*t120.*(t9.*t13.*v_e.*1.152921504606847e+18-t9.*t100.*v_e.*4.627631787495115e+31+t6.*t13.*v_n.*1.152921504606847e+18-t6.*t99.*v_n.*4.690217797021143e+31+h.*t9.*t103.*v_e.*1.460862238981533e+25+h.*t6.*t102.*v_n.*1.470707715282058e+25),-t37.*t59.*t108.*t120.*(et47+et48),t118];
mt6 = [0.0,0.0,t118.*v_d,t81+t4.*t118.*v_e,t118.*v_n.*-2.0,(t9.*t37.*t39.*t91)./(h.*t39.*t91+3.921472212741887e+33),t6.*t96.*t123,t9.*t10.*t37.*t96.*t123,0.0,t32.*t107,0.0,-t81-t4.*t118.*v_e.*2.0,t107.*v_d+t4.*t118.*v_n,-t80-t109.*2.0];
mt7 = [-t37.*t136.*(h.*t6.*(t18.*4.41939415613864e+33+t20.*1.489224625905806e+31+3.278727125438313e+35)+t6.*t39.*t90.*1.960736106370943e+33)];
mt8 = [(h.*t9.*t82.*2.0-t9.*t39.*t90.*6.848516479516672e+15)./(t13.*t82.*2.0-4.690217797021143e+31),-(t3.*t118)./t67-(t10.*t68)./(h.*t69+t69.*t105),0.0,0.0,-1.0,t118.*v_n,t80+t109,0.0,0.0,0.0,0.0,0.0,0.0,0.0,a_y.*t71+a_z.*t73,-a_y.*t74-a_z.*t70,t7.*(a_y.*t5-a_z.*t8),t10.*t37.*t65,-t8.*w_y-t5.*w_z,t37.*t65,0.0,0.0,0.0,t6.*t75,t9.*t75,-a_x.*t7-a_y.*t8.*t10-a_z.*t5.*t10,t38.*t136.*(et31+et32+et33+et34+et35),0.0];
mt9 = [t38.*t136.*(et49+et50+et51+et52+et53),0.0,0.0,0.0,-a_y.*t70+a_z.*t74-a_x.*t7.*t9,-a_y.*t73+a_z.*t71+a_x.*t6.*t7,0.0];
mt10 = [(t37.*(et36+et37+et38+et39))./(t13.*3.32306998946229e+35+t111+t122+t13.*t52.*4.449178648656757e+33+t13.*t52.^2.*1.489224625905806e+31+1.342813353646476e+49)];
mt11 = [t67.*t68.*7.29211651452621e-5-t68.*t109-t9.*t118.*v_n,-t5.*t37.*(t74.*t125-t71.*t118.*v_n)+t8.*t37.*(t70.*t125-t73.*t118.*v_n)];
out1 = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11],9,9);
end

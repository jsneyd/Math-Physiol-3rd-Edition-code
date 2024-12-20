
#fetal circulation
;
restart;
eq1:=Psa-Ps-Qs*Rsa;
eq2:=Ps-Psv-Qs*Rsv;
eq3:=Ppa-Ppv-Qp*Rp;
eq4:=F*(Cld*Ppv-Cls*Psa)-Ql;
eq5:=F*(Crd*Psv-Crs*Ppa)-Qr;
eq6:=Ql+Qd-Qs;
eq7:=Qd+Qp-Qr;
eq8:=Qr+Qf-Qs;
eq9:=Csv*(Psv+Ps)/2-Vsv;
eq10:=Cp*(Ppa+Ppv)/2+V0p-Vp;
eq11:=Csa*(Psa+Ps)/2+V0s-Vsa;
eq13:=Ppa-Psa-Qd*Rd;
eq12:=Vsa+Vsv+Vp-Vt;
NULL;
Psv:=Ppv;
#simplifying;
Crs:=0;
Cls:=0;
solve({eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13},{Psa,Ppa,Ppv,Ps,Vsa,Vp,Vsv,Qr,Qf,Qs,Qd,Qp,Ql});
assign(%);
simplify(Qf/Qr);
simplify(Qd/Qr);
simplify(Qs/Qr);
simplify(Qp/Qr);
simplify(Ql/Qr);
NULL;
Crd:=Cd;
Cld:=Cd;
simplify(Qf/Qr);
simplify(Qd/Qr);
simplify(Qs/Qr);
simplify(Qp/Qr);
simplify(Ql/Qr);
NULL;

restart;
eq1:=Psa-Ps-Qs*Rsa;
eq2:=Ps-Psv-Qs*Rsv;
eq3:=Ppa-Ppv-Qp*Rp;
eq4:=F*(Cld*Ppv-Cls*Psa)-Ql;
eq5:=F*(Crd*Psv-Crs*Ppa)-Qr;
eq6:=Ql+Qd-Qs;
eq7:=Qd+Qp-Qr;
eq8:=Qr+Qf-Qs;
eq9:=Csv*(Psv+Ps)/2-Vsv;
eq10:=Cp*(Ppa+Ppv)/2+V0p-Vp;
eq11:=Csa*(Psa+Ps)/2+V0s-Vsa;
eq13:=Ppa-Psa-Qd*Rd;
eq12:=Vsa+Vsv+Vp-Vt;
Qf:=0;
solve({eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13},{Psa,Ppa,Ppv,Ps,Psv,Vsa,Vp,Vsv,Qr,Qs,Qd,Qp,Ql});
assign(%);

Crs:=0;
Cls:=0;
simplify(Ppv/Ql);
simplify(Psv/Qr)
;
factor(Psa/Ppa);
simplify(Ppv/Psv)
;
NULL;

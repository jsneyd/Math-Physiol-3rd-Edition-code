
restart;
eq1:=F*Cd*Pv-Q;
eq2:=Pa-Pv-Q*Rs;
eq3:=Ca*Pa-Va;
eq4:=Cv*Pv-Vv;
eq5:=Va+Vv-Vt;
solve({eq1,eq2,eq3,eq4,eq5},{Pa,Pv,Q,Va,Vv});
NULL;


restart;
km1 := k1*K1;
km2 := k2*K2;
km3 := k3*K3;
km4 := k4*K4;

eq1 := k2*pe - km2*pi + km3*si*ci - k3*pii;
eq2 := km2*pii - k2*pe + k1*se*ce - km1*pe;
eq3 := km4*ce - k4*ci + k3*pii - km3*si*ci;
eq4 := k4*ci - km4*ce + km1*pe - k1*se*ce;
eq5 := pii+pe+ci+ce-C0;

 solve([eq2,eq3,eq4,eq5],[ci,ce,pii,pe]);


assign(%);
J = simplify(-ci*km3*si + k3*pii);
NULL;

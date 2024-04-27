
restart;
eq1 := km1*Ki^2*x2 + k8*ATP*z1 - (km8+k1)*x1;        
eq2 := km2*x3 + k1*x1 - (km1*Ki^2 + k2*Ni^3)*x2;      
eq3 := km3*z2*ADP + k2*Ni^3*x2 - (k3+km2)*x3;         
eq4 := km4*y3 + k3*x3 - (k4 + km3*ADP)*z2;            
eq5 := k4*z2 + km5*Ne^3*y2 - (k5 + km4)*y3;          
eq6 := k5*y3 +km6*y1 - (km5*Ne^3+k6*Ke^2)*y2;         
eq7 := k6*Ke^2*y2 + km7*Pi*z1 - (km6 + k7)*y1;       
eq8 := x1 + x2 + x3 + y1 + y2 + y3 + z1 + z2 - 1;
solve([eq1 ,eq2 ,eq3 ,eq4 ,eq5 ,eq6 ,eq7 ,eq8 ],[x1,x2,x3,y1,y2,y3,z1,z2]);
assign(%);
J := k4*z2 - km4*y3;
simplify(J)
;
km1 := K1*k1;
km2 := K2*k2;
km3 := K3*k3;
km4 := K4*k4;
km5 := K5*k5;
km6 := K6*k6;
km7 := K7*k7;
km8 := K8*k8;
simplify(numer(J));


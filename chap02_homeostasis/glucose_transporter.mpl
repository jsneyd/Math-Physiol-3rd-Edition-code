
restart;
eq1 := k*pe - k*p_i + kp*si*ci - km*p_i;
eq2 := k*p_i - k*pe + kp*se*ce - km*pe;
eq3 := k*ce - k*ci + km*p_i - kp*si*ci;
eq4 := k*ci - k*ce + km*pe - kp*se*ce;
eq5 := p_i+pe+ci+ce-C0;
 solve({eq2,eq3,eq4,eq5},{ci,ce,p_i,pe})
;
assign(%)
;

J := simplify(km*p_i - kp*si*ci);
NULL;

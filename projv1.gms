$ondollar
$Title Projectv1
$offsymxref offsymlist offuelxref offuellist offupper
option limrow = 0, limcol = 0;

set rows /r1*r12/;
set cols /c1*c12/;
alias(rows, i,j,k,l);

parameter DSM(i,j);
parameter size;

size=4

 ;

DSM(i,j) = uniform(0,1);
DSM(i,i) = 0;


parameter BDSM(i,j)/
 r1.r3 1
 r3.r2 1
 r4.(r5,r6,r12) 1
 r5.(r6,r8,r11) 1
 r6.(r2,r9,r12) 1
 r7.(r2,r11) 1
 r8.(r1,r4,r9,r11) 1
 r9.(r3,r6,r10) 1
 r10.(r2,r3,r11,r12) 1
 r11.(r2,r3) 1
 r12.(r1,r9,r10,r11) 1

/;



*BDSM(i,j)$(DSM(i,j) >= 0.5) = 1;
*BDSM(i,j)$(DSM(i,j) < 0.5) = 0;
display BDSM;

free variables obj;

binary variable x(i), y(i,j);

equations objective1,
          cons1;

cons1..
sum(i, x(i)) =e= size;

objective1..
obj =e= sum((i,j), BDSM(i,j) * x(i) * x(j)) - 0* sum((i,j), BDSM(i,j) * x(i) * (1 - x(j))) + 0* sum((i,j), BDSM(i,j) * (1 - x(i)) * x(j));



model project1 /cons1, objective1/;


solve project1 using minlp maximizing obj;
display x.l;

parameter index(i);
parameter count;
count = 1;

loop(i,
if(x.l(i) = 1,
index(j)$(ord(j) = count) = ord(i);
count = count + 1;
)
);

loop(i,
if(x.l(i) = 0,
index(j)$(ord(j) = count) = ord(i);
count = count + 1;
)
);

parameter newBDSM(i,j);

loop(k,
loop(l,
newBDSM(i,j)$(ord(k) = index(i) and ord(l) = index(j)) = BDSM(k,l);
)
);
display index, BDSM, newBDSM;

equations
          objective2,
          cons2(j),
          cons3(i);
objective2..
obj =e= sum((i,j)$(ord(i) <=size and ord(j) <= size) , newBDSM(i,j) * y(i,j)) ;

cons2(j)$(ord(j) <= size)..
 sum(i$(ord(i) <= size), y(i,j)) =e=1;

cons3(i)$(ord(i)<=size)..
 sum(j$(ord(j) <= size), y(i,j)) =e=1;

 model project1a /objective2, cons2, cons3/;


solve project1a using minlp maximizing obj;


file soln /out.txt/;

*y=w*x+gamma

put soln;

loop(i,
loop(j,
put BDSM(i,j):2:0;
)
put /;
);

put /;
put /;

loop(i,
loop(j,
put newBDSM(i,j):2:0;
)
put /;
);

putclose soln;

execute 'gdxxrw @out.txt'

display obj.l, x.l, y.l;

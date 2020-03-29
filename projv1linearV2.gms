$ondollar
$Title Projectv1
$offsymxref offsymlist offuelxref offuellist offupper
option limrow = 0, limcol = 0;

set rows /r1*r12/;
set cols /c1*c12/;
alias(rows, i,j,k,l);

parameter DSM(i,j);
parameter size;

size=1

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

parameter BDSM1(i,j);
BDSM1(i,j) = BDSM(i,j);

*BDSM(i,j)$(DSM(i,j) >= 0.5) = 1;
*BDSM(i,j)$(DSM(i,j) < 0.5) = 0;
display BDSM;
parameter nextIndex "starting position of the next optimal block";
nextIndex = 1;

free variables obj;

binary variable x(i), y(i,j),z(i,j), u(i,j);

equations objective1,
          cons1
          cons2
          cons3(i,j)
          cons4(i,j)
          cons5
          cons6(i,j)
          cons7(i,j)
          cons8
          cons9(i,j)
          cons10(i,j);

cons1..
sum(i, x(i)) =e= size;

cons2..
sum((i,j), y(i,j)) =e= size * size;

cons3(i,j)..
y(i,j) =l= x(i);

cons4(i,j)..
y(i,j) =l= x(j);

cons5..
sum((i,j), z(i,j)) =e= size * (card(i) - size);

cons6(i,j)..
z(i,j) =l= x(i);

cons7(i,j)..
z(i,j) =l= (1-x(j));

cons8..
sum((i,j), u(i,j)) =e= size * (card(i) - size);

cons9(i,j)..
u(i,j) =l= (1-x(i));

cons10(i,j)..
u(i,j) =l= x(j);


objective1..
obj =e= sum((i,j), BDSM1(i,j) * y(i,j)) -
sum((i,j)$(ord(j) >= nextIndex and ord(i) >= nextIndex), BDSM1(i,j) * z(i,j))+
0*sum((i,j)$(ord(j) >= nextIndex and ord(i) >= nextIndex), BDSM1(i,j) * u(i,j));



model project1 /cons1, cons2, cons3, cons4,cons5, cons6, cons7, cons8, cons9, cons10, objective1/;
parameter index(i);
parameter count;
parameter newBDSM(i,j);

*loop(rows$(ord(rows) <= 2),
loop(rows,
x.fx(i)$(ord(i)<nextIndex) = 0;

solve project1 using mip maximizing obj;


*display x.l;


count = 1;

loop(i,
if(ord(i) < nextIndex,
index(i) = ord(i);
count = count + 1;
)
if(x.l(i) = 1,
index(j)$(ord(j) = count) = ord(i);
count = count + 1;
)
);

loop(i,
if(x.l(i) = 0 and ord(i) >= nextIndex,
index(j)$(ord(j) = count) = ord(i);
count = count + 1;
)
);




loop(k,
loop(l,
newBDSM(i,j)$(ord(k) = index(i) and ord(l) = index(j)) = BDSM1(k,l);
)
);

BDSM1(i,j) = newBDSM(i,j);
nextIndex = nextIndex + 1;
);



display index, BDSM, newBDSM;
$ontext
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
$offtext

file soln /outlin.txt/;

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

execute 'gdxxrw @outlin.txt'

display obj.l, x.l, y.l;

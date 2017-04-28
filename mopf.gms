*** 1.Sets definitions

Sets k periods of time    / ini,1*3/
     i  loop iterations   /1*10/
     n  buses             /b1,b2/
     t  thermal units     /a1,a2/
     nc(t,n)   bus units  / a1.b1,a2.b2/
     cgt thermal data     /pmax,pmin,vc,fc,sc/
     cln line data        /ss,cp/
     p(k) dynamic subset of time periods
     it(i) dynamic subset od iterations ;

alias (n,n1);  alias (k,j);

*** 2. Data presentation

Table G(t,cgt) thermal unit data

         pmax     pmin    vc         fc     sc
*        (puW)    (puW)   (pu$/puW)  (pu$)  (pu$)
a1       1.3      0.02    0.1        10.    20.
a2       2.5      0.02    0.125      17.    18.

Table  Ld(n,k) Load

        1      2       3
*       (puW)  (puW)   (puW)
b1      1.0    2.6     0.2


Table  Ln(n,n1,cln) line suscepptance and capacity limit
          ss         cp
*         (puOhm)    (puMVa)

b1.b2     1.2        1.5;


Scalars
          sp   upper bound cost   /inf/
          in lower bound cost     /-inf/;

Parameters
           l(i,t,k) marginal cost      / 1.(a1,a2).1*3  0./
           val(i,t,k) value of status variable / 1.(a1,a2).1*3 0./
           tot(i) subproblem cost   /1 0./;

*** 3.Declaration of variables
Positive Variables o(t,k) ,v(t,k) ,y(t,k),r;
Free Variables ag(n,k) ,objs,objm;
Binary Variables tpv(t,k) ,tpy(t,k) ;


*** 4. Declaration and definition of subproblem equations
Equations Cost         Benders subproblem cost
          Bpw(n,k)     Power balance at every bus
          Pma(t,k)     Maximum Output Power
          Pmi(t,k)     Minimum output Power
          Cma(n,n1,k)  n-n1 mac line transmission capacity
          Cmi(n,n1,k)  n1-n max line trapsmission capacity;


Cost(p)..objs=e=sum(t,G(t,"vc")*o(t,p));
Bpw(n,p)..sum(t$nc(t,n),o(t,p))+sum(n1,(Ln(n,n1,'ss')*(ag(n1,p)-ag(n,p))+Ln(n1,n,"ss")*(ag(n1,p)-ag(n,p))))=e=Ld(n,p);

Pma(t,p).. o(t,p)=l=G(t,'pmax')*v(t,p);
Pmi(t,p).. o(t,p)=g=G(t,'pmin')*v(t,p);
Cma(n,n1,p)..Ln(n,n1,'ss')*(ag(n,p)-ag(n1,p))=l=Ln(n,n1,'cp');
Cmi(n,n1,p)..Ln(n,n1,'ss')*(ag(n,p)-ag(n1,p))=g=Ln(n,n1,'cp');


*** 5. Declaration and definition of master problem equations
Equations Mast Benders master problem cost
Cut(i) Benders cut
Rvy(t,k) Logic of plant operation
Fact(k) Feasibility cut;


Mast..     objm=e=r+sum((k,t)$p(k),G(t,'fc')*tpv(t,k)+G(t,'sc')*tpy(t,k));
Cut(it)..  r=g=tot(it)+sum((t,k)$p(k),l(it,t,k)*(tpv(t,k)-val(it ,t,k))) ;
Rvy(t,k)$p(k)..tpy(t,k)=g=tpv(t,k)-tpv(t,k-1);

Fact(k)$p(k)..sum(t,G(t ,'pmax')*tpv(t,k))=g=sum(n,Ld(n,k));


*** 6 . Initial values
tpv.fx(t,'ini')=0;


*** 7. Model and Solve statement in a LOOP
Model SBP Benders Subproblem/Cost,Bpw,Pma,Pmi,Cma,Cmi/;
Model MST Benders Master Problem /Mast,Cut,Rvy,Fact/;


Loop(i$((abs(sp-in)/(abs(in)+1)) gt 1e-10),

* initial period is not considered
p(k) =yes;
p(k)$(ord(k) eq 1) = no;


* add a Benders cut and solve master problem
it(i) = yes;
Solve MST using MIP minimizing objm;


* variables updating
in=objm.L;
v.fx(t,k) = tpv.l(t,k);
y.fx(t,k) = tpy.L(t,k);
tot(i+1) = 0.;


* only first period is considered
p(k)$(ord(k) ne 2)= no;

  Loop(j$(ord(j)  gt 1),
* Solve subproblem J
  Solve SBP using LP minimizing objs;


* variables updating
  tot(i+1)     =  tot(i+1)+ objs.L;
  l(i+1,t,j)   =  v.m(t,j);
  val(i+1,t,j) =  v.l(t,j);



* next period is considered
 p(j)   =  no;
 p(j+1) =  yes;

 );

  sp=tot(i+1)+objm.l-r.l;

);

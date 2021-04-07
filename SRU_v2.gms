option mip = cplex;
option optcr  = 1e-5,  
       optca  = 0.5,
	 limrow = 500,
	 limcol = 500,
       reslim = 80000;
$include data.txt
*threads to use for computation
*mode of mipstart
*solve the dual first
*try the mipstart file 3 times with repair before moving on
*$onecho > cplex.opt
*threads 0
*mipstart 5
*predual 1
*repairtries 3
*$offecho

alias(n,nn);
alias(n,i);
alias(n,j);
alias(u,uu);
Binary variable y(n);
positive variable z(n,nn);
variable obj;
equation dummy; dummy   ..  obj=e=sum((n,nn)$(ord(n) gt ord(nn)),z(n,nn));
*next equation use parameter to select number of nodes equal to max_Num
equation total_nodes; total_nodes .. sum(n,y(n)) =e= max_Num;
*equation list_uniqueness(u); list_uniqueness(u) .. sum(n$list(u,n),y(n))=e=1;
*$include weighted_eqn.txt
equation list_weights(u); list_weights(u) .. sum(n$list(u,n),y(n))=e=weight(u);
equation connectivitya(n,nn); connectivitya(n,nn)$(ord(n) ne ord(nn)) .. z(n,nn) =l= y(n)*M(n,nn); 
equation connectivityb(n,nn); connectivityb(n,nn)$(ord(n) ne ord(nn)) .. z(n,nn) =g= M(n,nn)*(y(n)+y(nn)-1); 
equation connectivityas(n,nn); connectivityas(n,nn)$(ord(n) ne ord(nn)) .. z(n,nn) =l= y(nn)*M(n,nn); 

$include init_guess.txt
Model nonsharp / all /;
nonsharp.OptFile = 1;
solve nonsharp using mip maximizing obj;
execute_unload 'SRU_v2' y.l;
execute "gdx2sqlite -i SRU_v2.gdx -o SRU_v2.db";


*y.lo(n)=0;
*h.lo(i,l)=0;
*fl.lo(u)=0;
*ll.lo(u)=0;
*fl_ll.lo(u,u)=0;
*z.lo(n,nn)=0;

*y.up(n)=1;
*h.up(i,l)=1;
*fl.up(u)=1;
*ll.up(u)=1;
*fl_ll.up(u,u)=1;
*z.up(n,nn)=INF;

*nonsharp.OptFile = 1;
*solve nonsharp using mip maximizing obj;
*execute_unload 'SRU_v2' y.l;
*execute "gdx2sqlite -i SRU_v2.gdx -o SRU_v2.db";
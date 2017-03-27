a=md();
assert(a.val==[] && a.var==[]);
a=md(1); a=md(0); a=md(-5);
assert(a.var==0 && a.std==0);
b=md(a);
assert(b.val==a.val && b.var==a.var);
assert(isa(b,'md'))

%%
syms x y z
f=symfun(x*exp(y)+z,[x y z]);
a=md([1 2 3],0.1);
md.exprInc(f,a)
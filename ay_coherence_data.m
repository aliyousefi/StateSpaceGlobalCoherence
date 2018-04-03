rho_a = 0.7;rho_b = 0.6;rho_c = 0.8;

%%
nd = 32;
U  = orth(randn(nd,nd));
la = diag(exp(-rho_a*(0:nd-1)));   
Sa = U*la*U';
la = svd(Sa);
ga=la(1)/sum(la);

lb = diag(exp(-rho_b*(0:nd-1)));
Sb = U*lb*U';
lb = svd(Sb);
gb=lb(1)/sum(lb);

lc = diag(exp(-rho_c*(0:nd-1)));
Sc = U*lc*U';
lc = svd(Sc);
gc=lc(1)/sum(lc);
%% generate three 3 sessions of data
df = 32;
Gs = [];
for i=1:100000
    S = wishrnd(Sa,df);
    ls = svd(S);
    Gs = [Gs;ls(1)/sum(ls)];
end
for i=1:100
    S = wishrnd(Sb,df);
    ls = svd(S);
    Gs = [Gs;ls(1)/sum(ls)];
end
for i=1:100
    S = wishrnd(Sc,df);
    ls = svd(S);
    Gs = [Gs;ls(1)/sum(ls)];
end

true_g = [ga*ones(100,1);gb*ones(100,1);gc*ones(100,1)];
plot(Gs);hold on;plot(true_g);hold off


%%
cl;
lambda = 1000;
alpha = 2.5;
phi = 300;

t = lambda^(-alpha);

x = 1:2000;



y = alpha*t*(x-phi).^(alpha-1).*exp(-t*(x-phi).^alpha);

myplot(1:2000, y);


%% 



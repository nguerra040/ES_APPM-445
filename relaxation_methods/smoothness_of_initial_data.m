% Smoothness of data
N = 16;
x = 1:N;
h = 1/(N+1);

% Data (a)
u = (-1).^x;
figure(1)
plot(u)

% Data (b)
u = sin(8*pi.*x*h);
figure(2)
plot(u)

% Data (c)
u = sin(pi.*x*h);
figure(3)
plot(u)

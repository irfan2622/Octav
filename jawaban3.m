clc;% untuk membersihkan jendela perintah atau konsol
clear all;% menghapus semua variabel dari ruang kerja
close all;% menutup semua jendela gambar atau plot


function [x, t, u] = courant(d, f, lb, xb, xu, tb, tu, dx, dt) % Mendefinisikan fungsi courant(upwind)
  t = tb:dt:tu;
  x = xb:dx:xu;
  u = [];
  
  for j = 1:length(x)
    u(j, 1) = f(x(j));
  endfor 
  
  for n = 1:length(t)
    u(1, n) = lb(t(n));
  endfor 
  
  c = d * dt / dx;
  for n = 1:length(t)-1
    for j = 2:length(x)
      u(j, n+1) = (1-c) * u(j, n) + c * u(j-1, n);
    endfor
  endfor
endfunction

function [x, t, u] = richardson(d, f, lb, rb, xb, xu, tb, tu, dx, dt) % Mendefinisikan fungsi richardson
  t = tb:dt:tu;
  x = xb:dx:xu;
  nt = length(t);
  nx = length(x);
  u = [];
  
  for j = 1:nx
    u(j, 1) = f(x(j));
  endfor
  
  for n = 1:nt
    u(1, n) = lb(t(n));
    u(nx, n) = rb(t(n));
  endfor
  
  c = (d*dt) / (2*dx);
  for n = 1:nt-1
    for j = 2:nx-1
      u(j, n+1) = u(j, n) - (c * (u(j+1, n) - u(j-1, n)));
    endfor
  endfor
endfunction

function [x, t, u] = lax(d, f, lb, rb, xb, xu, tb, tu, dx, dt) % Mendefinisikan fungsi lax
  t = tb:dt:tu;
  x = xb:dx:xu;
  nt = length(t);
  nx = length(x);
  u = [];
  
  for j = 1:nx
    u(j, 1) = f(x(j));
  endfor
  
  for n = 1:nt
    u(1, n) = lb(t(n));
    u(nx, n) = rb(t(n));
  endfor
  
  c = (d*dt) / (2*dx);
  for n = 1:nt-1
    for j = 2:nx-1
      u(j, n+1) = (u(j+1, n) + u(j-1, n))/2 - (c * (u(j+1, n) - u(j-1, n)));
    endfor
  endfor
endfunction

d = 2; # d adalah nilai 𝑑 pada 𝑢𝑡 + 𝑑𝑢𝑥 = 0.
f = @(x) e^(-x.^2); # f adalah nilai awal f(x) pada persamaan transport
lb = @(t) e^(-(2*t).^2); # syarat batas kiri pada persamaan transport
rb = @(t) e^(-(2-2*t).^2) # syarat batas kanan pada persamaan transport
xb = 0; # batas bawah untuk variabel x
xu = 2; # batas atas untuk variabel x
tb = 0; # batas bawah untuk variabel t
tu = 1; # batas atas untuk variabel t
dx = 0.5; # stepsize x
dt = 0.1; # stepsize t

# hasil aproksimasi
[x, t, u1] = courant(d, f, lb, xb, xu, tb, tu, dx, dt);
[x, t, u2] = richardson(d, f, lb, rb, xb, xu, tb, tu, dx, dt);
[x, t, u3] = lax(d, f, lb, rb, xb, xu, tb, tu, dx, dt);

# print
[x'] # nilai x
[t'] # nilai t
[u1'] # hasil aproksimasi metode courant
[u2'] # hasil aproksimasi metode richardson
[u3'] # hasil aproksimasi metode lax

sol = @(x, t) e^(-1*(x-2*t)^2); # solusi eksak
for j = 1:length(x)
  for n = 1:length(t)
    y(j, n) = sol(x(j), t(n));
  endfor
endfor

# Plot 1 (Metode Upwind)
figure(1);
mesh(x, t, u1');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Metode Upwind");

# Plot 2 (Metode Richardson)
figure(2);
mesh(x, t, u2');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Metode Richardson");

# Plot 3 (Metode Lax)
figure(3);
mesh(x, t, u3');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Metode Lax");

# Plot 4 (Solusi Eksak)
figure(4);
mesh(x, t, y');
xlabel('x');
ylabel('t');
zlabel('u');
title("Solusi Eksak")


clc; % untuk membersihkan jendela perintah atau konsol
clear all; % menghapus semua variabel dari ruang kerja
close all; % menutup semua jendela gambar atau plot

function [x_i, w_1i, w_2i] = nonlinshoot(f, fy, fyp, a, b, n, alpha, beta, m, tol)  % Mendefinisikan fungsi nonlinearshoot 
  h = (b - a)/n;% stepsize
  k = 1; % parameter k
  tk = (beta - alpha)/(b - a); % parameter tk
  x_i = w_1i = w_2i = []; % mendefinisikan variabel variabel xi,w1i dan w2i sebagai sebuah himpunan kosong
  while k <= m % while condition agar fungsi berjalan
    w = [alpha;tk]; % Mendefinisikan w(1,0) dan w(2,0)
    u = [0,1]; % Mendefinisikan nilai u1 dan u2
    for i = 1:n 
      x = a + (i-1)*h; % Mendefinisikan nilai x
      % proses rungkuta orde 4
      k_11 = h*w(2,i);
      k_12 = h*f(x, w(1,i), w(2,i));

      k_21 = h*(w(2,i)+(k_12/2));
      k_22 = h*f((x+(h/2)), (w(1,i)+(k_11/2)), (w(2,i)+(k_12/2)));

      k_31 = h*(w(2,i)+(k_22/2));
      k_32 = h*f((x+(h/2)), (w(1,i)+(k_21/2)), (w(2,i)+(k_22/2)));

      k_41 = h*(w(2,i)+k_32);
      k_42 = h*f((x+h), (w(1,i)+k_31), (w(2,i)+k_32));

      w(1,i+1) = w(1,i) + ((k_11 + 2*k_21 + 2*k_31 + k_41)/6);
      w(2,i+1) = w(2,i) + ((k_12 + 2*k_22 + 2*k_32 + k_42)/6);

      kp_11 = h*u(2);
      kp_12 = h*(fy(x, w(1,i), w(2,i))*u(1) + fyp(x, w(1,i), w(2,i))*u(2));

      kp_21 = h*(u(2) + (kp_12/2));
      kp_22 = h*(fy((x+(h/2)), w(1,i), w(2,i))*u(1) + fyp((x+(h/2)), w(1,i), w(2,i))*(u(2) + (kp_12/2)));

      kp_31 = h*(u(2)+(kp_22/2));
      kp_32 = h*(fy((x+(h/2)), w(1,i), w(2,i))*(u(1) + (kp_21/2)) + fyp((x+(h/2)), w(1,i), w(2,i))*(u(2) + (kp_22/2)));

      kp_41 = h*(u(2)+kp_32);
      kp_42 = h*(fy((x+h), w(1,i), w(2,i))*(u(1)+kp_31) + fyp((x+h), w(1,i), w(2,i))*(u(2) + kp_32));

      u(1) = u(1) + (kp_11 + 2*kp_21 + 2*kp_31 + kp_41)/6;
      u(2) = u(2) + (kp_12 + 2*kp_22 + 2*kp_32 + kp_42)/6;
    endfor

  if abs(w(1,n+1) - beta) <= tol       % jika sudah mencapai batas toleransi maka program berhenti
    for i = 1:(n+1)
      x = a+(i-1)*h;
      x_i(i) = x;
      w_1i(i) = w(1,i);
      w_2i(i) = w(2,i);
    endfor
    return
  endif
  tk = tk-((w(1,n+1) - beta)/u(1));
  k = k + 1;
  endwhile
  disp('max iteration')
endfunction

xi = w1i = w2i = []; %mendefinisikan variabel variabel xi,w1i dan w2i sebagai sebuah himpunan kosong
f = @(x, y, yp) ((1/2)*(1-((yp) ^ 2)- y*sin(x))); %mendefinisikan fungsi
fy = @(x, y, yp) (-(sin(x)) / 2); % melakukan turunan terhadap y
fyp = @(x, y, yp) (-yp); % melakukan turunan terhadap yp
a = 0; % mendefinisikan batas bawah fungsi
b = pi; % mendefinisikan batas atas fungsi
alpha = 2; % Mendefinisikan nilai awal fungsi
beta = 2; % Mendefinisikan nilai batas fungsi
tol = 10^(-4); % batas toleransi error
n = 20; % mendefinisikan parameter n
m = 10; % mendefinisikan parameter m

[xi, w1i, w2i] = nonlinshoot(f, fy, fyp, a, b, n, alpha, beta, m, tol); % Mendefinisikan fungsi nonlinearshoot 

% Mencari solusi eksak
sln = @(x) (2 + sin(x));
w = [];
for i = 1:length(xi)
  w(i) = sln(xi(i));
endfor

[xi', w1i', w', abs(w-w1i)'] % Mencetak tabelx,aproksimasi w, solusi eksak dan eror


hold on; 
fplot(sln, [0,pi], 'k'); % plot fungsi eksaknya
scatter(xi, w1i, 'r'); % xi-wli scatter plot  
legend('Solusi Eksak', 'Metode Nonlinear Shooting'); % memberi legenda pada plot
title('Nonlinear shooting method'); % judul plot


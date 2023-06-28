clc;% untuk membersihkan jendela perintah atau konsol
clear all;% menghapus semua variabel dari ruang kerja
close all;% menutup semua jendela gambar atau plot

function [t_grid,w]=nonlinear_FDM_naive(f,f_y,f_yp,a,b,n,alpha,beta,max_iter,TOL) % Mendefinisikan beda hingga nonlinear
  h=(b-a)/(n+1); %stepsize
  w=zeros(n,1); %vektor solusi aproksimasi
  t_grid=[a:h:b]; %Membuat mesh_point
  J=zeros(n,n); %Membuat matriks jacobian
  F=zeros(n,1); %Membuat vektor fungsi  F=(f_1,f_2,...,f_n) yang dievaluasi di x_k

  for i=1:n %inisialisasi solusi awal
    w(i)=alpha+i*(beta-alpha)/(b-a)*h;
  endfor
  k=1;
  while k<=max_iter %lakukan iterasi jika masih belum didapat kriteria stopnya

    %solve nonlinear sistem tersebut dengan metode newton
    x=a+h;
    %kontruksi matriks Jacobian, dan vektor F-nya
    t=(w(2)-alpha)/(2*h);
    J(1,1)=2+h^2*f_y(x,w(1),t); %main diagoanal
    J(1,2)=-1+(h/2)*f_yp(x,w(1),t); %right diagonal
    F(1)=(2*w(1)-w(2)-alpha+h^2*f(x,w(1),t));
    for i =2:n-1
      x=a+i*h;
      t=(w(i+1)-w(i-1))/(2*h);
      J(i,i)=2+h^2*f_y(x,w(i),t); %main diagoanal
      J(i,i+1)=-1+(h/2)*f_yp(x,w(i),t); %main diagoanal
      J(i,i-1)=-1-(h/2)*f_yp(x,w(i),t); %left diagoanal
      F(i)=(2*w(i)-w(i+1)-w(i-1)+h^2*f(x,w(i),t));
    endfor
     x=b-h;
     t=(beta-w(n-1))/(2*h);
     J(n,n)=2+h^2*f_y(x,w(n),t); %main diagonal
     J(n,n-1)=-1-(h/2)*f_yp(x,w(n),t); %right diagonal
     F(n)=(2*w(n)-w(n-1)-beta+h^2*f(x,w(n),t));



    v=inverse(J)*F; %vector v adalah product dari J^-1 F
    w= w-v; % lakukan update nilai pada w

    if norm(v,2)<= TOL %kriteria stop jika norm(v)<=toleransinya
      break;
     else
        k=k+1; %jika belum memenuhi kriteria stop terus lanjut iterasinya (memperbaiki nilai w)
    endif
  endwhile
  w=[alpha ; w ; beta]; %konstruksi akhir w
  t_grid=transpose(t_grid); % %transpose meshpoint
  % untuk konsistensi dimensi saja

endfunction


f=@(x,y,yp) yp+2*(y-log(x))-1/x ; %mendefinisikan fungsi
f_y=@(x,y,yp) 2*(y-log(x)); %melakukan turunan terhadap y
f_yp=@(x,y,yp) 1; %melakukan turunan terhadap yp
a=2; %mendefinisikan batas bawah fungsi
b=3; %mendefinisikan batas atas fungsi
alpha=0.5+log(2); %Mendefinisikan nilai awal fungsi
beta=1/3+ log(3); %Mendefinisikan nilai batas fungsi
n=9; %banyaknya partisi (pilih n=9 sehingga h=0.1)
max_iter=30; %Memilih maksimal iterasi
TOL=10^(-4); %batas toleransi error

%memanggil fungsi nonlinear_FDM_naive
[x_grid,w]=nonlinear_FDM_naive(f,f_y,f_yp,a,b,n,alpha,beta,max_iter,TOL)
f_anal= @(x) 1./x +log(x); %sol analitik

%membuat grafiknya
fplot(f_anal, [a,b],'b')  % plot fungsi eksaknya
hold on;
scatter(x_grid,w,'r') % scatter plot  
legend('solusi analitik', 'solusi linear FDM'); % memberi legenda pada plot
legend("location", "northwest");



%membuat tabel saja.
sol_anal=f_anal(x_grid); %sol analitik di meshpoint
error=abs(w-sol_anal); %error
[x_grid,w,sol_anal,error]
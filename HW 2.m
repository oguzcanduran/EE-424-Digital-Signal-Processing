n=[0:255];
x_n=zeros(1,256);
x_n(1:256)=cos((pi/256).*n.^2);

figure();
stem(x_n);
title("x[n] vs n")
xlabel("n")
ylabel("x[n]")


h_0=zeros(1,256);
h_0(1)=1/sqrt(2);
h_0(2)=1/sqrt(2);

h_1=zeros(1,256);
h_1(1)=1/sqrt(2);
h_1(2)=-1/sqrt(2);

H0_z=fft(h_0);
H1_z=fft(h_1);
F0_z=flip(H0_z);
F1_z=-flip(H1_z);

X0_z=H0_z.*fft(x_n);
X1_z=H1_z.*fft(x_n);


f_1=ifft(F1_z);
f_0=ifft(F0_z);

x0n=ifft(X0_z);
x1n=ifft(X1_z);

%x0n=conv(h_0,x_n);
%x1n=conv(x_n,h_1);

figure();
stem(x0n);
title("x_0[n]")
xlabel("n")
ylabel("x_0[n]")


figure();
stem(x1n);
title("x_1[n]")
xlabel("n")
ylabel("x_1[n]")

v0=zeros(1,128);
for n =1:128
    v0(n)=x0n(2*n);
end

figure();
stem(v0);
title("v_0[n]")
xlabel("n")
ylabel("v_0[n]")

v1=zeros(1,128);
for n =1:128
    v1(n)=x1n(2*n);
end

figure();
stem(v1);
title("v_1[n]")
xlabel("n")
ylabel("v_1[n]")

y0=zeros(1,256);
for n =1:256
    if mod(n,2)==0
        y0(n)=v0(n./2);
    else
        y0(n)=0;
    end
end

figure();
stem(y0);
title("y_0[n]")
xlabel("n")
ylabel("y_0[n]")

y1=zeros(1,256);
for n =1:256
    if mod(n,2)==0
        y1(n)=v1(n/2);
    else
        y1(n)=0;
    end
end

figure();
stem(y1);
title("y_1[n]")
xlabel("n")
ylabel("y_1[n]")

Y0_z=fft(y0(1:256));
Y1_z=fft(y1(1:256));
sum=Y0_z.*F0_z+ Y1_z.*F1_z;
y=ifft(sum);

figure();
stem(y);
title("y[n] vs n")
xlabel("n")
ylabel("y[n]")





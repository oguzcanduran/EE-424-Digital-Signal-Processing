close all
rng(424); % get same results at every run
%Generate xn
n_1=0:15;
xn= normrnd(0,1,[1,16])+1j*normrnd(0,1,[1,16]);
disp("x[n]:")
disp(xn)

%Plotting x[n]
figure();
subplot(1,2,1);
stem([0:15],real(xn));
grid on;
xlabel("n");
ylabel("Re \{x[n]\}");
title("x[n] Real Part");

subplot(1,2,2);
stem([0:15],imag(xn));
grid on;
xlabel("n");
ylabel("Im \{x[n]\}");
title("x[n] Imaginary Part");

%% DFT of x[n]
XK=discrete_ft(xn);
disp("X(K)")
disp(XK)

figure();
subplot(1,2,1);
stem([0:15],real(XK));
grid on;
xlabel("K");
ylabel("Re \{X(K)\}");
title("X(K) Real Part");

subplot(1,2,2);
stem([0:15],imag(XK));
grid on;
xlabel("K");
ylabel("Im \{X(K)\}");
title("X(K) Imaginary Part");

%% Part a
x1n = xn(16:- 1:1) ;
figure();
subplot(1,2,1);
stem([0:15],real(x1n));
grid on;
xlabel("n");
ylabel("Re \{x_1[n]\}");
title("x_1[n] Real Part");

subplot(1,2,2);
stem([0:15],imag(x1n));
grid on;
xlabel("n");
ylabel("Im \{x_1[n]\}");
title("x_1[n] Imaginary Part");

X1K_1= discrete_ft(x1n); % Direct DFT
X1K_2= [XK(1,1) XK(16:- 1:2)].*exp(1j*(2*pi/16).*n_1); % Numerically
disp("Direct DFT :");
disp(X1K_1);
disp("Analytical Result: ");
disp(X1K_2);

figure();
subplot(2,2,1);
stem([0:15],real(X1K_1));
grid on;
xlabel("K");
ylabel("Re \{X_1(K)\}");
title("Directly Calculated X_1(K) Real Part");

subplot(2,2,2);
stem([0:15],imag(X1K_1));
grid on;
xlabel("K");
ylabel("Im \{X_1(K)\}");
title("Directly Calculated X_1(K) Imaginary Part");


subplot(2,2,3);
stem([0:15],real(X1K_2));
grid on;
xlabel("K");
ylabel("Re \{X_1(K)\}");
title("Analytically Calculated X_1(K) Real Part");

subplot(2,2,4);
stem([0:15],imag(X1K_2));
grid on;
xlabel("K");
ylabel("Im \{X_1(K)\}");
title("Analytically Calculated X_1(K) Imaginary Part");


%% Part b
x2n=((-1).^n_1).*xn;
figure();
subplot(1,2,1);
stem([0:15],real(x2n));
grid on;
xlabel("n");
ylabel("Re \{x_2[n]\}");
title("x_2[n] Real Part");

subplot(1,2,2);
stem([0:15],imag(x2n));
grid on;
xlabel("n");
ylabel("Im \{x_2[n]\}");
title("x_2[n] Imaginary Part");

X2K_1= discrete_ft(x2n); % Direct DFT
nb=mod(n_1+length(xn)/2,length(xn))+1;
X2K_2= XK([nb]); % Numerically
disp("Direct DFT :");
disp(X2K_1);
disp("Analytical Result: ");
disp(X2K_2);

figure();
subplot(2,2,1);
stem([0:15],real(X2K_1));
grid on;
xlabel("K");
ylabel("Re \{X_2(K)\}");
title("Directly Calculated X_2(K) Real Part");

subplot(2,2,2);
stem([0:15],imag(X2K_1));
grid on;
xlabel("K");
ylabel("Im \{X_2(K)\}");
title("Directly Calculated X_2(K) Imaginary Part");


subplot(2,2,3);
stem([0:15],real(X2K_2));
grid on;
xlabel("K");
ylabel("Re \{X_2(K)\}");
title("Analytically Calculated X_2(K) Real Part");

subplot(2,2,4);
stem([0:15],imag(X2K_2));
grid on;
xlabel("K");
ylabel("Im \{X_2(K)\}");
title("Analytically Calculated X_2(K) Imaginary Part");


%% Part c
% it is seen that two x[n] is added back to back and x3[n] is formed 
x3n=[xn xn];

figure();
subplot(1,2,1);
stem([0:31],real(x3n));
grid on;
xlabel("n");
ylabel("Re \{x_3[n]\}");
title("x_3[n] Real Part");

subplot(1,2,2);
stem([0:31],imag(x3n));
grid on;
xlabel("n");
ylabel("Im \{x_3[n]\}");
title("x_3[n] Imaginary Part");

X3K_1= discrete_ft(x3n); % Direct DFT
X3K_2=zeros(1,32); % Numerically
X3K_2(1:2:31)=2*XK;
disp("Direct DFT :");
disp(X3K_1);
disp("Analytical Result: ");
disp(X3K_2);


figure();
subplot(2,2,1);
stem([0:31],real(X3K_1));
grid on;
xlabel("K");
ylabel("Re \{X_3(K)\}");
title("Directly Calculated X_3(K) Real Part");

subplot(2,2,2);
stem([0:31],imag(X3K_1));
grid on;
xlabel("K");
ylabel("Im \{X_3(K)\}");
title("Directly Calculated X_3(K) Imaginary Part");


subplot(2,2,3);
stem([0:31],real(X3K_2));
grid on;
xlabel("K");
ylabel("Re \{X_3(K)\}");
title("Analytically Calculated X_3(K) Real Part");

subplot(2,2,4);
stem([0:31],imag(X3K_2));
grid on;
xlabel("K");
ylabel("Im \{X_3(K)\}");
title("Analytically Calculated X_3(K) Imaginary Part");

%% Part d

x4n=xn(1:8)+xn(9:16);

figure();
subplot(1,2,1);
stem([0:7],real(x4n));
grid on;
xlabel("n");
ylabel("Re \{x_4[n]\}");
title("x_4[n] Real Part");

subplot(1,2,2);
stem([0:7],imag(x4n));
grid on;
xlabel("n");
ylabel("Im \{x_4[n]\}");
title("x_4[n] Imaginary Part");

X4K_1= discrete_ft(x4n); % Direct DFT
n_half=0:7;
X4K_2=XK(n_half.*2+1);% Numerically
disp("Direct DFT :");
disp(X4K_1);
disp("Analytical Result: ");
disp(X4K_2);


figure();
subplot(2,2,1);
stem([0:7],real(X4K_1));
grid on;
xlabel("K");
ylabel("Re \{X_4(K)\}");
title("Directly Calculated X_4(K) Real Part");

subplot(2,2,2);
stem([0:7],imag(X4K_1));
grid on;
xlabel("K");
ylabel("Im \{X_4(K)\}");
title("Directly Calculated X_4(K) Imaginary Part");


subplot(2,2,3);
stem([0:7],real(X4K_2));
grid on;
xlabel("K");
ylabel("Re \{X_4(K)\}");
title("Analytically Calculated X_4(K) Real Part");

subplot(2,2,4);
stem([0:7],imag(X4K_2));
grid on;
xlabel("K");
ylabel("Im \{X_4(K)\}");
title("Analytically Calculated X_4(K) Imaginary Part");


%% Part e
x5n=zeros(1,32);
x5n(1:16)=xn;

figure();
subplot(1,2,1);
stem([0:31],real(x5n));
grid on;
xlabel("n");
ylabel("Re \{x_5[n]\}");
title("x_5[n] Real Part");

subplot(1,2,2);
stem([0:31],imag(x5n));
grid on;
xlabel("n");
ylabel("Im \{x_5[n]\}");
title("x_5[n] Imaginary Part");

X5K_1= discrete_ft(x5n); % Direct DFT
temp=discrete_ft(xn.*exp((-1j*pi/16).*n_1));
X5K_2=zeros(1,32);
X5K_2(1:2:31)=XK;
X5K_2(2:2:32)=temp;
disp("Direct DFT :");
disp(X5K_1);
disp("Analytical Result: ");
disp(X5K_2);


figure();
subplot(2,2,1);
stem([0:31],real(X5K_1));
grid on;
xlabel("K");
ylabel("Re \{X_5(K)\}");
title("Directly Calculated X_5(K) Real Part");

subplot(2,2,2);
stem([0:31],imag(X5K_1));
grid on;
xlabel("K");
ylabel("Im \{X_5(K)\}");
title("Directly Calculated X_5(K) Imaginary Part");


subplot(2,2,3);
stem([0:31],real(X5K_2));
grid on;
xlabel("K");
ylabel("Re \{X_5(K)\}");
title("Analytically Calculated X_5(K) Real Part");

subplot(2,2,4);
stem([0:31],imag(X5K_2));
grid on;
xlabel("K");
ylabel("Im \{X_5(K)\}");
title("Analytically Calculated X_5(K) Imaginary Part");



%% Part f

x6n=zeros(1,32);
x6n(1:2:31)=xn;

figure();
subplot(1,2,1);
stem([0:31],real(x6n));
grid on;
xlabel("n");
ylabel("Re \{x_6[n]\}");
title("x_6[n] Real Part");

subplot(1,2,2);
stem([0:31],imag(x6n));
grid on;
xlabel("n");
ylabel("Im \{x_6[n]\}");
title("x_6[n] Imaginary Part");

X6K_1= discrete_ft(x6n); % Direct DFT
X6K_2=[XK XK];%Numerically
disp("Direct DFT :");
disp(X6K_1);
disp("Analytical Result: ");
disp(X6K_2);



figure();
subplot(2,2,1);
stem([0:31],real(X6K_1));
grid on;
xlabel("K");
ylabel("Re \{X_6(K)\}");
title("Directly Calculated X_6(K) Real Part");

subplot(2,2,2);
stem([0:31],imag(X6K_1));
grid on;
xlabel("K");
ylabel("Im \{X_6(K)\}");
title("Directly Calculated X_6(K) Imaginary Part");


subplot(2,2,3);
stem([0:31],real(X6K_2));
grid on;
xlabel("K");
ylabel("Re \{X_6(K)\}");
title("Analytically Calculated X_6(K) Real Part");

subplot(2,2,4);
stem([0:31],imag(X6K_2));
grid on;
xlabel("K");
ylabel("Im \{X_6(K)\}");
title("Analytically Calculated X_6(K) Imaginary Part");




%% Part g

x7n=xn(1:2:15);

figure();
subplot(1,2,1);
stem([0:7],real(x7n));
grid on;
xlabel("n");
ylabel("Re \{x_7[n]\}");
title("x_7[n] Real Part");

subplot(1,2,2);
stem([0:7],imag(x7n));
grid on;
xlabel("n");
ylabel("Im \{x_7[n]\}");
title("x_7[n] Imaginary Part");

X7K_1= discrete_ft(x7n); % Direct DFT
n_half=1:8;
X7K_2=(XK(n_half) + XK(n_half+8))/2  ;% Numerically
disp("Direct DFT :");
disp(X7K_1);
disp("Analytical Result: ");
disp(X7K_2);


figure();
subplot(2,2,1);
stem([0:7],real(X7K_1));
grid on;
xlabel("K");
ylabel("Re \{X_7(K)\}");
title("Directly Calculated X_7(K) Real Part");

subplot(2,2,2);
stem([0:7],imag(X7K_1));
grid on;
xlabel("K");
ylabel("Im \{X_7(K)\}");
title("Directly Calculated X_7(K) Imaginary Part");


subplot(2,2,3);
stem([0:7],real(X7K_2));
grid on;
xlabel("K");
ylabel("Re \{X_7(K)\}");
title("Analytically Calculated X_7(K) Real Part");

subplot(2,2,4);
stem([0:7],imag(X7K_2));
grid on;
xlabel("K");
ylabel("Im \{X_7(K)\}");
title("Analytically Calculated X_7(K) Imaginary Part");


%% DFT FUNCTION
function xk=discrete_ft(x)
    N=length(x);
    s=0:N-1;
    exps =exp(-1j*(2*pi/N)*s);
    xk=x*(exps.'.^s);
end



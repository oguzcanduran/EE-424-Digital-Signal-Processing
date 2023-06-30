%% BSF

z1=zero_create(0.99*exp(1j*pi*0.42));
z2=zero_create(0.99*exp(1j*pi*0.5));
z3=zero_create(0.99*exp(1j*pi*0.52));
z4=zero_create(0.99*exp(1j*pi*0.56));
zeros=[z1 z2 z3 z4];

z_zero=z-zeros;
h_n=poly(zeros);
disp(h_n)
n=0:length(h_n)-1;
figure();
zplane(h_n, 1);
ylabel('jIm\{z\}');
xlabel('Re\{z\}');
title("Pole-Zero plot of Band Stop Filter")

figure();
stem(n,h_n)
title('Impulse Response of the BSF')
ylabel('h[n]')
xlabel('n')

eq = 1;
for i = 1:length(z_zero)
eq= eq * z_zero(i);
end

w=-pi:pi*0.002:pi;
Hz=  sum(h_n.'.*exp(n.'.*(-1j).*w));
figure()
plot(omega, abs(Hz)/max(abs(Hz)));
title('Magnitude of Response of H_{BSF}(e^{jw})')
ylabel('∣ H_{BSF}(e^{jw} ∣')
xlabel('\omega')

figure();
freqz(sym2poly(eq)/max(abs(Hz)))

%% 3
k=0:1:1023;
x_f=cos(k.^2*(pi/512));
figure();
plot(k,x_f);
axis tight
xlabel("n")
ylabel("x_f[n]")
title("Chirp Signal Plot")
y_f=conv(x_f,h_n/max(abs(Hz)));
figure();
plot(y_f)
axis tight
xlabel("n")
ylabel("y_f[n]")
title("Chirp Signal After BSF Plot")

m = audioread("music.mp3");
filtered=conv(m(:,1),h_n/max(abs(Hz)));
audiowrite('music_bsf.wav', real(filtered), 44100)

record=audioread("ses.m4a");
filtered_2=conv(record,h_n/max(abs(Hz)));
audiowrite("record_bsf.wav", real(filtered_2), 48000)












function q=zero_create(z)
q=zeros(1,4);
q(1)=z;
q(2)=1/z;
q(3)=conj(z);
q(4)=conj(1/z);
end







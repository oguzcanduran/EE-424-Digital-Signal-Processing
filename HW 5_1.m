
clear all
close all
syms z

%% LPF
zerolar=zeros(4,39);
zerolar(:,1)=zero_create(0.99*exp(1j*pi*0.106));
for pp=1:39
zeros=zerolar(:,1:pp);
len=4*length(zeros(1,:));
z_zero=z-zeros;
h_n=real(poly([-1 reshape(zeros,[1,len])]));
n=0:length(h_n)-1;
eq = 1;
for i = 1:length(z_zero)
eq= eq * z_zero(i);
end

w=-pi:pi*0.002:pi;
Hz=  sum(h_n.'.*exp(n.'.*(-1j).*w));
plot(w,abs(Hz)/max(abs(Hz)))
temp=Hz(1:420);
ma=find(abs(temp)==max(abs(temp)));
disp(500-ma);
in=((500-ma)/500)*pi;
zerolar(:,pp+1)=zero_create(0.99*exp(1j*in));
end
figure();
plot(w,abs(Hz)/max(abs(Hz)))

function q=zero_create(z)
q=zeros(1,4);
q(1)=z;
q(2)=1/z;
q(3)=conj(z);
q(4)=conj(1/z);
end





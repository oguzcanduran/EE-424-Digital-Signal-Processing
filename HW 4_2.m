%% E
d=0:119;
hn=0.96.^d;
lxn=137216;
n=0:lxn-1;
xn=cos(0.05*pi.*n);
%xn=xn(:,1)';
l_xn=length(xn);
l_hn=length(hn);

y=zeros(1,l_xn+l_hn-1);
for q= 1:l_xn+l_hn-1

    if q <= l_hn-1
        y(q)= (exp(1j*pi*(q-1)/20)*(1-(0.96*exp(-1j*pi/20))^(q))/ ...
            (1-(0.96*exp(-1j*pi/20))) + exp(-1j*pi*(q-1)/20)* ...
            (1-(0.96*exp(1j*pi/20))^(q))/(1-(0.96*exp(1j*pi/20))))/2;
    elseif q>l_xn-1
        y(q)= (exp(1j*pi*(q-1)/20)*((0.96*exp(-1j*pi/20))^(q-l_xn)- ...
            (0.96*exp(-1j*pi/20))^(l_hn))/(1-(0.96*exp(-1j*pi/20))) + ...
            exp(-1j*pi*(q-1)/20)*((0.96*exp(1j*pi/20))^(q-l_xn)- ...
            (0.96*exp(1j*pi/20))^(l_hn))/(1-(0.96*exp(1j*pi/20))))/2;
    else 
        y(q)= (exp(1j*pi*(q-1)/20)*(1-(0.96*exp(-1j*pi/20))^(l_hn))/ ...
            (1-(0.96*exp(-1j*pi/20))) + ...
            exp(-1j*pi*(q-1)/20)*(1-(0.96*exp(1j*pi/20))^(l_hn))/ ...
            (1-(0.96*exp(1j*pi/20))))/2;
    end
end
y4=y(121:length(y));
%% EA
[y1,comp]=convolution(xn,hn);
y1=y1(121:length(y1));
%% EB
Ns=1024-l_hn+1; %Ns=905
xn_new=[xn zeros(1,Ns- mod(l_xn,Ns))];  % 344 zero added
l_xn_2=length(xn_new);
n_part=l_xn_2/Ns;
[comph, hk ]=dft(hn, 1024);
tot_comp=comph;
y2=zeros(1,l_xn_2+l_hn-1); 
for i=1:n_part
    [temp_comp1, XK]=dft(xn_new((i-1)*Ns+1:i*Ns),1024);
    [temp_comp2, yk]=idft(hk.*XK,1024);
    tot_comp=tot_comp+ temp_comp1 + temp_comp2 ;
    y2((i-1)*Ns+1:(i-1)*Ns+1024)= y2((i-1)*Ns+1:(i-1)*Ns+1024)+ yk;
end
y2=y2(121:length(y2)-344); 

%% EC
Nsc=1024-l_hn+1; %Ns=905
xn_newc=[xn zeros(1,Ns- mod(l_xn,Ns))];  % 344 zero added
l_xn_2c=length(xn_newc);
n_partc=l_xn_2c/Nsc;
hkc=fft(hn, 1024);
y3=zeros(1,l_xn_2+l_hn-1);
for i=1:n_partc
    XKc=fft(xn_newc((i-1)*Ns+1:i*Ns),1024);
    ykc=ifft(hkc.*XKc,1024);
    y3((i-1)*Ns+1:(i-1)*Ns+1024)= y3((i-1)*Ns+1:(i-1)*Ns+1024)+ ykc;
end
y3=y3(121:length(y3)-344);
%% ED
[loss14,e14]=loss(y1,y4);
[loss24,e24]=loss(y2,y4);
[loss34,e34]=loss(y3,y4);

%% EF
figure();
stem((20000:20511),xn(20000:20511));
title('Plot of x[n] for 19999<n<20512 )');
xlabel('n');
ylabel('x[n]');
axis tight
grid on
figure();
stem((20000:20511),y1(20000:20511));
title('Plot of y_1[n] for 19999<n<20511 )');
xlabel('n');
ylabel('y_1[n]');
axis tight
grid on
%% FUNCTIONS
% PART A FUNCTION
function [y1,num]=convolution(x,h)
lx=length(x);
lh=length(h);
y1=zeros(1,lx+lh-1);
num=0;
for n= 2: lx+lh
    
    y=0;
    for k=1:lh
        if n-k<=0
            xtemp=0 ;
        elseif n-k >lx
            xtemp=0;
        else
            xtemp=x(n-k);
        end
        y=y+ h(k)*xtemp;
        num=num+1;
    end
    y1(n-1)=y;
end
end 
% PART B FUNCTION
function [complexity, xk]=dft(x,N) % N indicates the N point dft
    if length(x)<N
        z=zeros(1,N);
        z(1:length(x))=x;
        x=z;
    end
    s=0:N-1;
    exps =exp(-1j*(2*pi/N)*s);
    xk=x(1:N)*(exps.'.^s);  % multiplications are done in this stage
    % this stage above is a matrix multiplication stage
    % the dimensions are (1,N) and (N,N)
    % so a row from left is multiplied with a column from left
    % so for each column N multiplications are needed
    % since there are N different columns
    % there are N^2 multiplications
    complexity=N^2;

end

function [complexity, xk]=idft(x,N) % N indicates the N point dft
    if length(x)<N
        z=zeros(1,N);
        z(1:length(x))=x;
        x=z;
    end
    s=0:N-1;
    exps =exp(1j*(2*pi/N)*s);
    xk= (1/N)*(x(1:N)*(exps.'.^s));
    %similar to dft the complexity is N^2
    complexity=N^2;
end
% PART D FUNCTION
function [loss,max_ep]=loss(y1,y2)
    N=length(y1);
    loss=0;
    max_ep=0;
    for i= 1:N
        dist=y1(i)-y2(i);
        loss=loss + dist^2/N;
        if max_ep<abs(dist)
            max_ep=dist;
        end
    end
end


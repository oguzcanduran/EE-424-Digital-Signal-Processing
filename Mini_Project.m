close all
clear all

%% PRELIMINARIES

k(1,:)=[1,0,1,5,4];
k(2,:)=[0,1,1,3,8];
N=256;
n=0:N-1;
exps_p=exp(1j*(2*pi/N).*n);

row=exps_p'.^k(1,:);
col=exps_p'.^k(2,:);
for i=1:5
    for j=1:256
        Y(i,j,:)=row(:,i).*col(j,i);
    end
    y_disp=round(line_eq(real(Y(i,:,:))));
    y_disp=reshape(y_disp,[N,N]);
    figure();
    imshow(uint8(y_disp));
    title(sprintf('k_1=%d ,k_2=%d', ...
        k(1,i),k(2,i)))
    imwrite(uint8(y_disp),sprintf([ ...
        '1-a-i-%d.tif'], i));
end
%% SIZE REDUCTION
image=imread("ProjectPictureMod.tif");

temp_image=image;
down_results=zeros(3,4096,4096);
for j=1:3
    down=size_reduction(temp_image,4);
    temp_image=down;
    s=size(down);
    down_results(j,1:s(1),1:s(2))=down;
    figure();
    imshow(down);
    title(sprintf('Downsampled to %d x %d', ...
        4096/(4^j),4096/(4^j)))
    imwrite(down,sprintf(['2-a-i-%d.tif' ...
        ], j));
end

temp_image=image;
lpf_down_results=zeros(3,4096,4096);
dft_coef=zeros(size(image));
for j=1:3
    dft_image=dft_2d(temp_image);
    if j==1
        dft_coef=dft_image;   
       % for part 3 approximations)
    end

    lpf=lp_filter(temp_image,4);
    lpf_image=dft_image.*lpf;
    lpf_down=real(idft_2d(lpf_image));
    if j==1
        lpf_down=uint8(lpf_down);
        lpf_down=round(line_eq(lpf_down));
        figure();
        imshow(lpf_down);
        title('Filtered Image')
        imwrite(lpf_down,'2-b-i-1.tif');
    end

    lpf_down=uint8(size_reduction(lpf_down,4));
    lpf_down=round(line_eq(lpf_down));
    temp_image=lpf_down;
    s=size(lpf_down);
    lpf_down_results(j,1:s(1),1:s(2))=lpf_down;
    figure();
    imshow(lpf_down);
    title(sprintf(['Filtered and Downsampled ' ...
        'to %d x %d'], 4096/(4^j),4096/(4^j)))
    imwrite(lpf_down,sprintf([ ...
        '2-b-i-%d .tif'], j+1) );
end
%% APPROXIMATIONS

dft_image=dft_coef;
abs_image=abs(dft_image);
flat=reshape(abs(dft_image) ,4096*4096,1);
[temp,inds]=sort(flat,"descend");
magnitudes=temp(1:20);
indices=zeros(20,2);

for i=1:20
    
    row=mod(inds(i),4096);
    col=fix(inds(i)/4096)+1;
    if row==0
        col=col-1;
        row=4096;
    end
    indices(i,:)=[row-1 col-1];
end
disp(indices);

for j=1:4
    th_mat=threshold(dft_image,4^j);
    app_im=ifft2(th_mat,"symmetric");
    display_image=uint8(app_im);
    display_image=round(line_eq(display_image));
    figure();
    imshow(display_image);
    title(sprintf(['%d/%d Non-Zero DFT' ...
        ' Coefficients'],1, 4^j) )
    imwrite(display_image,sprintf(['' ...
        '3-b-i-%d.tif'],j))
end



%% FUNCTIONS
function y=line_eq(mat)
ma=max(max(mat));
mi=min(min(mat));

a=255/(ma-mi);
b=-a*mi;
y=a*mat+b;
end

function y=size_reduction(im,M)
    shape=size(im);
    r=shape(1);
    c=shape(2);
    y=im(1:M:r,1:M:c);
end



function lpf=lp_filter(x,r)
    len=length(x);
    len_reduce=len/(r*2);
    lpf=ones(len,len);
    lpf(len_reduce+1:len-len_reduce,:)=0;
    lpf(:,len_reduce+1:len-len_reduce)=0;
end

function y=dft_2d(mat)
shape=size(mat);
r=shape(1);
c=shape(2);
dft_rows=zeros(r,c);
y=zeros(r,c);
    for i=1:r
        dft_rows(i,:)=fft(mat(i,:));
    end
    for j=1:c
        y(:,j)=fft(dft_rows(:,j));
    end
end

function y=idft_2d(mat)
shape=size(mat);
r=shape(1);
c=shape(2);
idft_rows=zeros(r,c);
y=zeros(r,c);
    for i=1:r
        idft_rows(i,:)=ifft(mat(i,:));
    end
    for j=1:c
        y(:,j)=ifft(idft_rows(:,j));
    end
end

function TH=threshold(imag,M)
    abs_im=abs(imag);
    shape=size(imag);
    r=shape(1);
    c=shape(2);
    fla=reshape(abs_im,r*c,1);
    [sorted,inds]=sort(fla,"descend");
    t_len=length(sorted);
    q=reshape(imag,r*c,1);
    q(inds((t_len/M)+1:t_len))=0; 
    TH=reshape(q,r,c);
end

clc;
clear;
close all;
utilsPath = '../../../utils';
addpath(utilsPath);
addpath('../../../label');
tic
%% setup parameters 
H =10; % tuning forbetter performance
L=3; % tuning forbetter performance
C=10;%10
p=0.8;
k=5;
%% input data
% strDir='D:\博士资料\[开源整理]（孙杨师兄代码）\[MFD-NMoG]\MFDNMoG\NoisyDataLT3\15\';
strDir='D:\博士资料\[开源整理]\ASTTV-NTLA公开\data\';
for i=1:k
    picname=[strDir  num2str(i),'.bmp'];
    I=imread(picname);
    [~, ~, ch]=size(I);
    if ch==3
        I=rgb2gray(I);
    end
    D(:,:,i)=I;
end
        tenD=double(D);
        size_D=size(tenD);
        [n1,n2,n3]=size(tenD);
        T=C*sqrt(n1*n2);
        n_1=max(n1,n2);%n(1)
        n_2=min(n1,n2);%n(2)
        patch_frames=L;%L
        patch_num=n3/patch_frames;
for l=1:patch_num
    for i=1:patch_frames
        X(:,:,i)=tenD(:,:,patch_frames*(l-1)+i);
    end
%% initialization
        mu = 1e-2;
        lambda1=0.005;
        lambda2= H / sqrt((n_1*patch_frames));        
        lambda3=100;%100
        weight=[1,1,1];
        [tenB,tenT,tenN] = ATVSTIPT1LPLS(X, lambda1,lambda2,lambda3, weight,mu,T,p);
        for i=1:patch_frames
            tarImg=tenT(:,:,i);
            a=uint8(tarImg);
            figure,imshow(a,[]);
            backImg=tenB(:,:,i);b=uint8(backImg);
        end 
end
toc
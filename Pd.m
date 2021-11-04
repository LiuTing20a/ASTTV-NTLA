function [counter] = Pd(tarImg,location,m,n,l,tol)
%Pd 统计一幅图像目标检测的正确数目
%input
%tarImg 输入分割后的目标图像
%l 真实目标个数
%location 真实目标位置，需要导入目标位置矩阵，维度为l*2，l为真实目标个数，第一列为x坐标，第二列为y坐标
%m ,n 目标真实大小
% tol 分割目标阈值
% output counter 正确被检测的目标数目
[M,N]=size(tarImg);
%% 生成 0-1 label矩阵
label=zeros(M,N);
for i=1:M
    for j=1:N
        if tarImg(i,j)>=tol
            label(i,j)=1;
        end
    end
end
%% 提取目标位置
counter=0;
z=1;
temp=zeros(l,m*n);
for k=1:l
    for i=location(k,1):location(k,1)+m-1
        for j=location(k,2):location(k,2)+n-1
            temp(k,z)=label(i,j);
            z=z+1;
        end
    end
    z=1;
    k=k+1;
end
for i=1:l
    if max(temp(i,:))>0
        counter=counter+1;
    end
end


function [counter_fa] = Fa(tarImg,location,m,n,l,tol)
%Fa 统计一幅图像虚假目标检测的数目
%input
%tarImg 输入分割后的目标图像
%l 真实目标个数
%location 真实目标位置，需要导入目标位置矩阵，维度为l*2，l为真实目标个数，第一列为x坐标，第二列为y坐标
%m ,n 目标真实大小
% tol 分割目标阈值
% counter_fa 虚假目标数目
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
%% 让目标处直接为0
for k=1:l
    for i=location(k,1):location(k,1)+m-1
        for j=location(k,2):location(k,2)+n-1
            label(i,j)=0;
        end
    end
    k=k+1;
end
counter_fa=sum(sum(label~=0));
% if tol==0
%     counter_fa=0;
% end
% % 判断是否在3*3邻域内
% [row,col]=find(label);
% num=length(row);
% d=zeros(1,num*(num-1)/2);z=zeros(1,num*(num-1)/2);
% k=1;
% for i=1:num-1
%     for j=i+1:num
%     d(k)=abs(row(i)-row(j));
%     k=k+1;
%     end
% end
% k=1;
% for i=1:num-1
%     for j=i+1:num
%     z(k)=abs(col(i)-col(j));
%     k=k+1;
%     end
% end
% for i=1:length(d)
%     if d(i)<2 && z(i)<2
%         counter_fa=counter_fa-1;
%     end
% end



function [counter_fa] = Fa(tarImg,location,m,n,l,tol)
%Fa ͳ��һ��ͼ�����Ŀ�������Ŀ
%input
%tarImg ����ָ���Ŀ��ͼ��
%l ��ʵĿ�����
%location ��ʵĿ��λ�ã���Ҫ����Ŀ��λ�þ���ά��Ϊl*2��lΪ��ʵĿ���������һ��Ϊx���꣬�ڶ���Ϊy����
%m ,n Ŀ����ʵ��С
% tol �ָ�Ŀ����ֵ
% counter_fa ���Ŀ����Ŀ
[M,N]=size(tarImg);
%% ���� 0-1 label����
label=zeros(M,N);
for i=1:M
    for j=1:N
        if tarImg(i,j)>=tol
            label(i,j)=1;
        end
    end
end
%% ��Ŀ�괦ֱ��Ϊ0
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
% % �ж��Ƿ���3*3������
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



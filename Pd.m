function [counter] = Pd(tarImg,location,m,n,l,tol)
%Pd ͳ��һ��ͼ��Ŀ�������ȷ��Ŀ
%input
%tarImg ����ָ���Ŀ��ͼ��
%l ��ʵĿ�����
%location ��ʵĿ��λ�ã���Ҫ����Ŀ��λ�þ���ά��Ϊl*2��lΪ��ʵĿ���������һ��Ϊx���꣬�ڶ���Ϊy����
%m ,n Ŀ����ʵ��С
% tol �ָ�Ŀ����ֵ
% output counter ��ȷ������Ŀ����Ŀ
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
%% ��ȡĿ��λ��
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


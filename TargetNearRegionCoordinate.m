function [Coordinates_x, Coordinates_y]= TargetNearRegionCoordinate(x, y, Diameter)
%����Ŀ���ٽ���������
%Ŀ������[x,y]
%Ŀ�������СDiameter,����
Radius = (Diameter-1)/2;
Coordinates_x = [1:Diameter]-Radius-1 +x;
Coordinates_y = [1:Diameter]-Radius-1 +y;
end


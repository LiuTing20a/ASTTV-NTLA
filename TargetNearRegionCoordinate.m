function [Coordinates_x, Coordinates_y]= TargetNearRegionCoordinate(x, y, Diameter)
%计算目标临近区域坐标
%目标坐标[x,y]
%目标邻域大小Diameter,奇数
Radius = (Diameter-1)/2;
Coordinates_x = [1:Diameter]-Radius-1 +x;
Coordinates_y = [1:Diameter]-Radius-1 +y;
end


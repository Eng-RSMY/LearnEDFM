function [row,col] = intersectionsSegments(XY1,XY2)
%SEGMENT_INTERSECTIONS Intersections of line segments.

row = [];
col = [];

for i=1:length(XY1)
    for j=1:length(XY2)
        s1_x = XY1(i,3) - XY1(i,1);     s1_y = XY1(i,4) - XY1(i,2);
        s2_x = XY1(j,3) - XY1(j,1);     s2_y = XY1(j,4) - XY1(j,2);

        s = (-s1_y * (XY1(i,1) - XY1(j,1)) + s1_x * (XY1(i,2) - XY1(j,2))) / (-s2_x * s1_y + s1_x * s2_y);
        t = ( s2_x * (XY1(i,2) - XY1(j,2)) - s2_y * (XY1(i,1) - XY1(j,1))) / (-s2_x * s1_y + s1_x * s2_y);

        if (s > 0 && s < 1 && t > 0 && t < 1) %if (s >= 0 && s <= 1 && t >= 0 && t <= 1) /*The commented version included merely touching segments. */
            touch = true;
        else
            touch = false; % No collision
        end
        if (touch)
            row = [row i j];
            col = [col j i];
        end
    end
end
      
end
    
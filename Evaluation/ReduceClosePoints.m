function C = ReduceClosePoints(C ,l_boundry)
%https://softwareengineering.stackexchange.com/questions/225443/how-to-reduce-close-points
    x_size = abs(l_boundry(1,1)-l_boundry(1,2)); %220
    DELTA = x_size/10;
    T = [];
    for i=1:size(C,2)
        sum = C(:,i);
        n=1;
        for j=1:size(C,2)
            if i~=j     %not same point
                D = DistancePointToPoint(C(:,i),C(:,j));
                if D < DELTA
                    sum = sum + C(:,j);
                    n=n+1;
                end
            end
        end
        sum = sum./n; %new point -> save in T matrix
        T = [T sum];
    end
    C = T;
end
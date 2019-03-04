function [ t1,t2 ] = find_back_angles( x,xa,h,ha )
    
    m = (((xa-x)^2+(ha-h)^2)^(0.5))/2;
    theta = atan((xa-x)/(ha-h));
    phi = acos(m/15);
    
    t2 = -theta - phi;
    t1 = -theta + phi;
end


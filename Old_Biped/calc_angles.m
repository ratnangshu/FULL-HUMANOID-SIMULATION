%This function takes 5 inputs (xa, ha, x, h, l) and gives 3 outputs (q3, q4, q5)

 function [q3, q4, q5] = calc_angles(xa, ha, x, h, l)
    x = -x;
    q3 = -2*atan((4*(l*x - l*xa))/(h^2 - 2*h*ha + 2*l*h + ha^2 - 2*l*ha + x^2 - 2*x*xa + xa^2) + ((-(h^2 - 2*h*ha + ha^2 + x^2 - 2*x*xa + xa^2)*(h^2 - 2*h*ha + ha^2 - 4*l^2 + x^2 - 2*x*xa + xa^2))^(1/2) - 2*l*x + 2*l*xa)/(h^2 - 2*h*ha + 2*l*h + ha^2 - 2*l*ha + x^2 - 2*x*xa + xa^2));
    q3 = -q3;
    q4 = - 2*atan((4*(l*x - l*xa))/(h^2 - 2*h*ha + 2*l*h + ha^2 - 2*l*ha + x^2 - 2*x*xa + xa^2) + ((-(h^2 - 2*h*ha + ha^2 + x^2 - 2*x*xa + xa^2)*(h^2 - 2*h*ha + ha^2 - 4*l^2 + x^2 - 2*x*xa + xa^2))^(1/2) - 2*l*x + 2*l*xa)/(h^2 - 2*h*ha + 2*l*h + ha^2 - 2*l*ha + x^2 - 2*x*xa + xa^2)) - 2*atan(((-(h^2 - 2*h*ha + ha^2 + x^2 - 2*x*xa + xa^2)*(h^2 - 2*h*ha + ha^2 - 4*l^2 + x^2 - 2*x*xa + xa^2))^(1/2) - 2*l*x + 2*l*xa)/(h^2 - 2*h*ha + 2*l*h + ha^2 - 2*l*ha + x^2 - 2*x*xa + xa^2));
    q4 = -q4;
    q5 = -2*atan(((-(h^2 - 2*h*ha + ha^2 + x^2 - 2*x*xa + xa^2)*(h^2 - 2*h*ha + ha^2 - 4*l^2 + x^2 - 2*x*xa + xa^2))^(1/2) - 2*l*x + 2*l*xa)/(h^2 - 2*h*ha + 2*l*h + ha^2 - 2*l*ha + x^2 - 2*x*xa + xa^2));
 end

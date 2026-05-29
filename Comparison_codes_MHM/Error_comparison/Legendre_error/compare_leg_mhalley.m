function [x_value, y_value, x1_value, y1_value]=compare_leg_mhalley(n)
    % Error calculation of Legendre polynomial using FOM
    [x_double, y_double]=legendre_mhalley(n);
    [x_vpa, y_vpa]=legendre_mhalley_vpa(n);

    % Relative calculation
    rel_error = abs(1 - x_double ./ x_vpa);
    % x value
    x_value = x_double; 

   for i=1:length(x_value)
        if rel_error(i)>0
            y_value(i)=log(rel_error(i));
        else
            y_value(i)=-50;
        end
   end


       % Relative calculation
    rel_error1 = abs(1 - y_double ./ y_vpa);
    % x value
    x1_value = y_double; 

   for i=1:length(x1_value)
        if rel_error1(i)>0
            y1_value(i)=log(rel_error1(i));
        else
            y1_value(i)=-50;
        end
   end

end
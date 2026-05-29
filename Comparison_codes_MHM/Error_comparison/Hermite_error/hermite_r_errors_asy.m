function [x_value, y_value, x1_value, y1_value] = hermite_r_errors_asy(n)
    % Ensure the Symbolic Math Toolbox is available
    assert(license('test','Symbolic_Toolbox'), 'Symbolic Toolbox is required');

    % Double precision ASY
    [x_double, y_double] = hermite_asy(n);  

    % Extended precision ASY
    [x_vpa_dbl, y_vpa_dbl] = hermite_asy_ext(n);    


    % Ensure dimensions match
    assert(length(x_double) == length(x_vpa_dbl), 'Output lengths differ');
    assert(length(y_double) == length(y_vpa_dbl), 'Output lengths differ');

    % Compute relative errors
    x_double=x_double(floor(n/2)+mod(n,2)+1:n);
    x_vpa_dbl=x_vpa_dbl(floor(n/2)+mod(n,2)+1:n);
    rel_error = abs(1 - x_double ./ x_vpa_dbl);
    x_value = x_double;  % your x values
for i=1:length(x_value)
    if rel_error(i)>0
     y_value(i)=log(rel_error(i));
    else
      y_value(i)=-50;
    end
end
y_double=y_double(floor(n/2)+mod(n,2)+1:n);
% % y_double=flip(y_double);
y_vpa_dbl= y_vpa_dbl(floor(n/2)+mod(n,2)+1:n);
% % y_vpa_dbl=flip(y_vpa_dbl);
rel_error1 = abs(1 - y_double ./ y_vpa_dbl); 
x1_value = y_double;
for i=1:length(x1_value)
    if rel_error1(i)>0
     y1_value(i)=log(rel_error1(i));
    else
      y1_value(i)=-50;
    end
end

end

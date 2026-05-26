clear all
clc;
n_values = [1000000, 1050000, 1100000, 1200000, 1300000];
num_runs = 10;

time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method3 = zeros(size(n_values));

for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = cpu_time_fourth_test_leg(n, num_runs);
    time_method2(j) = cpu_time_third_leg(n, num_runs);
    time_method3(j) = cpu_time_halley_leg(n, num_runs);
   
end

%% plotting the results

s1=time_method1;
s2=time_method2;
s3=time_method3;
s1 = s1'; 
s2 = s2';
s3 = s3';


for j = 1:length(n_values)
    n = round(n_values(j));    
    [f, I1(j), ~] = fsweep_fourth_1(n);
    k1(j) = I1(j);
    [t, I2(j), ~] = fsweep_third(n);
    k2(j) = I2(j);
    et(j) = max(abs(t-f)./f);
    [h, I3(j), ~] = fsweep_halley_1(n);
    k3(j) = I3(j);
    eh(j) = max(abs(h-f)./f);

end  
% n1 = sum(I1);
% n2 = sum(I2);
% n3 = sum(I3);

%% plotting the results

t1=k1;
t2=k2;
t3=k3;

t1 = t1';
t2 = t2';
t3 = t3';
 e1 = et';
 e2 = eh';


%%
%% table output
format long g
n_values = n_values';

data = [n_values, t1, s1, t2, e1, s2, t3, e2, s3];

% ---- Print header ----
fprintf('%8s | %20s | %30s | %30s\n','n','FOM-H','TOM-H','MHM');
fprintf('%8s | %10s %10s | %10s %10s %10s | %10s %10s %10s\n',...
'', 'I.Count','A.Time', ...
'I.Count','Error','A.Time', ...
'I.Count','Error','A.Time');

fprintf('%s\n', repmat('-',1,130));

% ---- Print data ----
for i = 1:size(data,1)
    fprintf('%8d | %10.0f %10.4f | %10.0f %10.2e %10.4f | %10.0f %10.2e %10.4f\n',...
        data(i,1), ...
        data(i,2), data(i,3), ...
        data(i,4), data(i,5), data(i,6), ...
        data(i,7), data(i,8), data(i,9));
end

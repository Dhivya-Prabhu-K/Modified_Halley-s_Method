function avg_time = cpu_time_third_leg(n, num_runs)
times = zeros(1, num_runs);

for run = 1:num_runs
    tic;
    legendre_third(n);
    times(run) = toc;
end

avg_time = mean(times);
end
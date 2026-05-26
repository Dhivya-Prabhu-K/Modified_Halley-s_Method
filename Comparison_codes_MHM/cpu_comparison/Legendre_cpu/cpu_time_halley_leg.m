function avg_time = cpu_time_halley_leg(n, num_runs)
times = zeros(1, num_runs);

for run = 1:num_runs
    tic;
    legendre_mhalley(n);
    times(run) = toc;
end

avg_time = mean(times);
end
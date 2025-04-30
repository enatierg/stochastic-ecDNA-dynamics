%Last Modified: 12/12/2024

clear all;
close all;

% setting the parameters for simulation for questions 6 to 14

s=1;                            % proliferation rate for cells with ecDNA
M=1000;                         % total number of simulations
t_final=10;                     % the max simulation time
max_k=200;                      % maximum number of ecDNA copies to track
                                % (the number was chosen after running
                                % multiple simulations)
max_cells=100000;               % large number
t_p=0:0.1:t_final;              % time points for M(1)(t) and M(2)(t)
                                % estimation
t_points=length(t_p);           % same as above

% pre-allocating array for the simulation loop
N=zeros(M,max_k+1,t_points);    % tracks the N_k(t)


% setting up the loop for stochastic simulation
for m = 1:M
    % initialising arrays and parameters
    cell_pop=zeros(1,max_cells); % tracks the cell population over the loop
    cell_pop(1) = 1;     % initialising with one cell containing one ecDNA
    cell_count = 1;      % tracks the actual number of cells
    t=0;                 % setting starting time to 0
    t_idx=1;

    while t<t_final
        % setting up the proliferation rates
        rates=1+(cell_pop(1:cell_count)>0)*(s-1);
        total_rate=sum(rates);

        % establishing time updates using gillespies
        dt=exprnd(1/total_rate);
        t=t+dt;

        % calculating Nk(t) at the correct time intervals
        while t_idx<t_points && t>=t_p(t_idx+1)
            % counting the amount of cells with each ecDNA count
            count=histcounts(cell_pop(1:cell_count),0:(max_k+1));
            N(m,:,t_idx)=count;
            % updating the time step
            t_idx=t_idx+1;
        end

        % randomly choosing a parent cell to divide via weighted sampling
        chosen_cell=randsample(1:cell_count, 1,true,rates);
        parent_ecDNA=cell_pop(chosen_cell);

        if parent_ecDNA>0
            % per the division process, doubling the ecDNA count within
            % parent cell
            doubled_ecDNA=2*parent_ecDNA;
            % using binomial distribution to allocate how many of ecDNA
            % copies go into the first daughter cell
            daughter1_ecDNA=binornd(doubled_ecDNA, 1/2);
            % with the rest going into the second daughter cell
            daughter2_ecDNA=doubled_ecDNA-daughter1_ecDNA;
        else
            % if the parent cell does not have any ecDNA, it just splits
            daughter1_ecDNA=0;
            daughter2_ecDNA=0;
        end

        % updating the cell population
        cell_pop(chosen_cell)=daughter1_ecDNA;
        cell_pop(cell_count+1)=daughter2_ecDNA;

        % increasing the number that tracks the total amount of cells
        cell_count=cell_count+1;
    end

    % since the simulation ends early, calculating N_k(t_final)
    if t_idx==t_points
        count=histcounts(cell_pop(1:cell_count),0:(max_k+1));
        N(m,:,t_idx)=count;
    end
    % adding (optional) simulation tracker, to see whether the simulations
    % are running ahead
    disp(['the simulation ',num2str(m), ' out of ', num2str(M), ' has been completed.']);
end

%% (Question 6)
% calculating the proportions over the entire simulated data
proportions=N./sum(N,2);
% calculating the average proportions over the data with t=10 
avg_proportions = mean(proportions(:,:,t_points),1);
% and std errors
std_err = std(proportions(:,:,t_points))/sqrt(M);


% plotting the average cell proportions pk(t) over the number of ecDNA
% copies k
figure(1);
bar(0:(length(avg_proportions)-1),avg_proportions);
% looking only at the first 20 samples as instructed
axis ([0 20 0 1]);
% adjusting the sides for better visibility
xlim([-1 21]);
xlabel('number of ecDNA copies (k)');
ylabel('average cell proportions p_{k}(t) at t=10');
title('Distribution of the cell proportions with s=1');
grid on;
hold on;
% adding error bars
errorbar(0:(length(avg_proportions)-1), avg_proportions,std_err,'k.', 'LineStyle','none');
hold off;

%% (Question 7 & 8)
% initialising arrayts to store M1 and M2 for the for loop
m1 = zeros(M, t_points);    % == M(1)(t)
m2 = zeros(M, t_points);    % == M(2)(t)
k_values = 0:max_k;         % vector consisting of k values

% we calculate using basic formulas by running through all simulations
for m = 1:M
    for t = 1:t_points
        % removing unneccessary dimensions
        p_k_t=squeeze(proportions(m,:,t));
        % using formula M(1)(t)=sum(k*pk(t)) to calculate first moment
        m1(m,t)=sum(k_values.*p_k_t);
        % using formula M(2)(t)=sum(k^2*pk(t)) to calculate second moment
        m2(m,t)=sum((k_values.^2).*p_k_t);
    end
end

% calculating the mean values
m1_mean=mean(m1, 1);
m2_mean=mean(m2,1);

% adding in theoretical prediction from question 4 & 5
m1_theoretical=1;
m2_theoretical=1+t_p;

% plotting the average first moment M(1)(t) over time (t)
figure(2);
plot(t_p, m1_mean,'blue','Linewidth', 1.2);

% adding the theoretical prediction for comparison
hold on;
plot(t_p,m1_theoretical*ones(1,t_points),'--k', 'Linewidth',1.5);
xlabel ('time (sec)');
ylabel('first moment M^{(1)}(t)');
title('Evolution of the first moment M^{(1)}(t)');
legend('Simulated M^{(1)}(t)', 'Analytical M^{(1)}(t)');
grid on;

% plotting the average second moment M(2)(t) over time (t)
figure(3);
plot(t_p, m2_mean,'blue','Linewidth', 1.2);

% adding the theoretical prediction for comparison
hold on;
plot(t_p, m2_theoretical,'--k','Linewidth', 1.5);
xlabel ('time (sec)');
ylabel('second moment M^{(2)}(t)');
title('Evolution of the second moment M^{(2)}(t)');
legend('Simulated M^{(2)}(t)', 'Analytical M^{(2)}(t)');
grid on;

%% (Question 11)
% summing across all simulations for the case k=0
N_minus=sum(N(:,1,:));

% creating a for loop that calculates the same for k>0
N_plus=0;
for k=2:max_k
    N_plus=N_plus+sum(N(:,k,:));
end

% calculating N_minus/N_plus
N_div=squeeze(N_minus./N_plus);
% popping down the given theoretical estimate
N_div_estimate=t_p/2;

% plotting N_div against time to show that ~t/2
figure(4);
plot(t_p,N_div,'blue','Linewidth',1.2)
hold on;
% and adding the estimate to the plot
plot(t_p, N_div_estimate,'--k','Linewidth', 1.5);
xlabel ('time (sec)');
ylabel('N^-/N^+');
title('Evolution of N^-/N^+ over time');
legend('Simulated N^-/N^+', 'Predicted N^-/N^+');
grid on;

%% (Question 13)
% using eq. (8) to calculate f_minus from the simulation data
f_minus=squeeze(N_minus./(N_minus+N_plus));

% calculating theoretical prediction as determined in question 12
f_minus_theoretical=1-2./(t_p+2);

% plotting f_minus over time with the theoeretical prediction
figure(5);
plot(t_p,f_minus,'blue','Linewidth', 1.2)
hold on;
% adding the theoeretical prediction line to the plot
plot(t_p, f_minus_theoretical,'--k','Linewidth', 1.5);
xlabel ('time (sec)');
ylabel('f^-');
title('Evolution of f^- over time');
legend('Simulated f^-', 'Predicted f^-');
grid on;

%% (Question 14)
% as before, calculating f_plus using eq. (8)
f_plus=squeeze(N_plus./(N_minus+N_plus));

% adding in the theoretical prediction as derived in question 14
f_plus_theoretical=2./(t_p+2);

% plotting f_plus over time
figure(6);
plot(t_p,f_plus,'blue','Linewidth', 1.2)
hold on;
% adding the theoeretical prediction line to the plot
plot(t_p, f_plus_theoretical,'--k','Linewidth', 1.5);
hold on;
xlabel ('time (sec)');
ylabel('f^-');
title('Evolution of f^+ over time');
legend('Simulated f^+', 'Predicted f^+');
grid on;

% Defining different time steps.
h_values = [0.1, 0.05, 0.01];

% Initial condition.
y0 = 2;

% Defining ODE and exact solution dy/dx = 4*e^(0.9x)
f = @(x) 4*exp(0.9*x);
exact_solution = @(x) 5*exp(0.9*x) - 3;

% Create a new figure with two subplots to show values and errors between valeus and exact solution

figure;
subplot(2,1,1); % first subplot for the solutions.
hold on;
subplot(2,1,2); % second subplot for the errors.
hold on;

% For loop over time steps ( take each time steps to integrate.)
for j = 1:length(h_values)
    h = h_values(j); % current time step
    x = 0:h:10; % time span
    N = length(x); % number of steps

    % Initializing variable
    y = zeros(1, N);
    y(1) = y0;

    % Implementing Runge-Kutta method
    for i = 1:(N-1)
        % Calculating k values
        k1 = h*f(x(i));
        k2 = h*f(x(i) + 0.5*h);
        k3 = h*f(x(i) + 0.5*h);
        k4 = h*f(x(i) + h);
        % Updating next y for each time with for loop. 
        y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % Calculating error as a percentage of the exact solution
    exact = exact_solution(x);
    error = abs((y - exact)./exact) * 100;

    % Plotting solutions
    subplot(2,1,1);
    plot(x, y);

    % Ploting errors
    subplot(2,1,2);
    plot(x, error);
end

% Adding legends and labels
subplot(2,1,1);
legend('h = 0.1', 'h = 0.05', 'h = 0.01');
xlabel('x');
ylabel('y');
title('Solution of the ODE using Runge-Kutta method with different time steps');

subplot(2,1,2);
legend('h = 0.1', 'h = 0.05', 'h = 0.01');
xlabel('x');
ylabel('Error (%)');
title('Error of the numerical solution compared to the exact solution % (as a percentage)');
hold off;

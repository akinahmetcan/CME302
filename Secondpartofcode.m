% Defining different time steps
h_values = [0.1, 0.05, 0.01];

% Initial conditions.
y0 = [4; 6];

% Defining ODEs
f1 = @(y1, y2) -0.3*y1;
f2 = @(y1, y2) 4 - 0.2*y1 - 0.1*y2;

% Creating a new figure with two subplots
figure;
subplot(2,1,1); % first subplot for y1
hold on;
subplot(2,1,2); % second subplot for y2
hold on;

% For loop over time steps
for j = 1:length(h_values)
    h = h_values(j); % current time step
    x = 0:h:10; % time span
    N = length(x); % number of steps

    % Initializing variables
    y1 = zeros(1, N);
    y2 = zeros(1, N);
    y1(1) = y0(1);
    y2(1) = y0(2);

    % Implement Runge-Kutta method
    for i = 1:(N-1)
        % Calculatig k values for y1
        k1 = h*f1(y1(i), y2(i));
        k2 = h*f1(y1(i) + 0.5*k1, y2(i));
        k3 = h*f1(y1(i) + 0.5*k2, y2(i));
        k4 = h*f1(y1(i) + k3, y2(i));
        % Update y1
        y1(i+1) = y1(i) + (k1 + 2*k2 + 2*k3 + k4)/6;

        % Calculating k values for y2
        k1 = h*f2(y1(i), y2(i));
        k2 = h*f2(y1(i), y2(i) + 0.5*k1);
        k3 = h*f2(y1(i), y2(i) + 0.5*k2);
        k4 = h*f2(y1(i), y2(i) + k3);
        % Updating y2
        y2(i+1) = y2(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % Plotting y1
    subplot(2,1,1);
    plot(x, y1);

    % Plotting y2
    subplot(2,1,2);
    plot(x, y2);
end

% Adding legends and labels
subplot(2,1,1);
legend('h = 0.1', 'h = 0.05', 'h = 0.01');
xlabel('x');
ylabel('y1');
title('Solution of the first ODE using Runge-Kutta method with different time steps');

subplot(2,1,2);
legend('h = 0.1', 'h = 0.05', 'h = 0.01');
xlabel('x');
ylabel('y2');
title('Solution of the second ODE using Runge-Kutta method with different time steps');
hold off;

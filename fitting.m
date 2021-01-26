% DOA based localization
% Two nearby transmitters, distant from a vehicle-based reciever
% Goal: Cluster observations from each path, localize each transmitter
% Method: Update DOA clusters by fitting to inverse tangent curves

% @author {wwhoward}@vt.edu
% wireless@vt

clc;clear;close all
addpath('munkres', 'data');

%% Set constants

saveStr = "0_0625rmse0_04sep";
%saveStr = "test";
saveFlg = 0;


rx.initial_x = 25000;   % initial reciever x position in meters
rx.initial_y = -15000;  % initial reciever y position in meters
rx.speed_x = 0;         % reciever x speed in m/s
rx.speed_y = 50;        % reciever y speed in m/s

tx1.x = 0;
tx1.y = 0;
tx2.x = 5000;
tx2.y = 0;

rmse = 0.5*pi/180;  % RMSE of DOA estimates in radians (from deg)

duration = 10;          % Collection duration in minutes
interval = 0.1;        % How often to measure in seconds
t = linspace(0, 60*duration, 60*duration/interval);

%% Simulation

rx.pos_x = rx.initial_x + t*rx.speed_x;
rx.pos_y = rx.initial_y + t*rx.speed_y;

DOA1 = atan2(rx.pos_y - tx1.y, rx.pos_x - tx1.x) + rmse*randn(size(t));
DOA2 = atan2(rx.pos_y - tx2.y, rx.pos_x - tx2.x) + rmse*randn(size(t));

DIS1 = sqrt((rx.pos_y - tx1.y).^2 + (rx.pos_x - tx1.x).^2);
DIS2 = sqrt((rx.pos_y - tx2.y).^2 + (rx.pos_x - tx2.x).^2);


% Shuffle DOA observations together
DOA = zeros(1,2*length(t));
DIS = zeros(1,2*length(t));
for i=1:length(t)
    flip = round(rand);
    if flip==1
        DOA(2*i-1) = DOA1(i);
        DOA(2*i)   = DOA2(i);
    else
        DOA(2*i-1) = DOA2(i);
        DOA(2*i)   = DOA1(i);
    end
end

DOA1_hat = zeros(size(t));
DOA2_hat = zeros(size(t));
a = zeros(2, length(t));
b = zeros(2, length(t));
a(:,1) = 250;
b(:,1) = 500;
fun = @(x, t)atan((t-x(1))/x(2));
a_true =  ([tx1.y, tx2.y] - rx.pos_y(1))./rx.speed_y;
b_true = -([tx1.x, tx2.x] - rx.pos_x(1))./rx.speed_y;

tx1.true_curve = fun([a_true(1),b_true(1)], t);
tx2.true_curve = fun([a_true(2),b_true(2)], t);
for i=1:length(t)
    % Form clusters
    %labels = kmeans(DOA(1:2*i)'-atan((t(1:i)-mean(a(:,i)))/mean(b(:,i))), 2);
    weights = abs([DOA(2*i-1:2*i)-fun([a(1,i), b(1,i)], t(i)); DOA(2*i-1:2*i)-fun([a(2,i), b(2,i)], t(i))]);
    labels = munkres(weights);
    
    if labels(1,1) == 1
        DOA1_hat(i) = DOA(2*i-1);
        DOA2_hat(i) = DOA(2*i);
    else
        DOA2_hat(i) = DOA(2*i-1);
        DOA1_hat(i) = DOA(2*i);
    end
    
    % Estimate parameters based on clustered DOA's
    %track1_hat = lsqcurvefit(fun, [a(1,i), b(1,i)], t(1:i), DOA1_hat(1:i), [], [], optimset('Algorithm','levenberg-marquardt','Display','off')); 
    %track2_hat = lsqcurvefit(fun, [a(2,i), b(2,i)], t(1:i), DOA2_hat(1:i), [], [], optimset('Algorithm','levenberg-marquardt','Display','off'));
    track1_hat = lsqcurvefit(fun, [a(1,i), b(1,i)], t(1:i), DOA1_hat(1:i), [], [], optimset('Display','off')); 
    track2_hat = lsqcurvefit(fun, [a(2,i), b(2,i)], t(1:i), DOA2_hat(1:i), [], [], optimset('Display','off'));
    
%     track1_hat = lsqcurvefit(fun, [a(1,i), b(1,i)], t(1:i), DOA1(1:i), [], [], optimset('Display','off')); 
%     track2_hat = lsqcurvefit(fun, [a(2,i), b(2,i)], t(1:i), DOA2(1:i), [], [], optimset('Display','off'));
    
    % Estimate parameters based on accurately associated DOA's
%     track1     = lsqcurvefit(fun, [a(3,i), b(3,i)], t(1:i), DOA1(1:i), [], [], optimset('Algorithm','levenberg-marquardt','Display','off'));
%     track2     = lsqcurvefit(fun, [a(4,i), b(4,i)], t(1:i), DOA2(1:i), [], [], optimset('Algorithm','levenberg-marquardt','Display','off'));
     %track1     = lsqcurvefit(fun, [a(3,i), b(3,i)], t(1:i), DOA1(1:i), [], [], optimset('Display','off'));
     %track2     = lsqcurvefit(fun, [a(4,i), b(4,i)], t(1:i), DOA2(1:i), [], [], optimset('Display','off'));

    
    a(1,i+1) = track1_hat(1);
    b(1,i+1) = track1_hat(2);
    a(2,i+1) = track2_hat(1);
    b(2,i+1) = track2_hat(2);
%     a(3,i+1) = track1(1);
%     b(3,i+1) = track1(2);
%     a(4,i+1) = track2(1);
%     b(4,i+1) = track2(2);
    
    tx1.y_hat(i) = a(1,i)*rx.speed_y + rx.pos_y(1);
    tx2.y_hat(i) = a(2,i)*rx.speed_y + rx.pos_y(1);
    tx1.x_hat(i) = rx.pos_x(1) - b(1,i)*rx.speed_y;
    tx2.x_hat(i) = rx.pos_x(1) - b(2,i)*rx.speed_y;
    
%     tx1.est_delta(i) = sqrt((tx1.y_hat(i)-tx1.y_hat(i+1))^2 + (tx1.x_hat(i)-tx1.x_hat(i+1))^2);
%     tx1.est_delta(i) = sqrt((tx2.y_hat(i)-tx2.y_hat(i+1))^2 + (tx2.x_hat(i)-tx2.x_hat(i+1))^2);

    tx1.error(i) = sqrt((tx1.x_hat(i)-tx1.x).^2 + (tx1.y_hat(i)-tx1.y).^2);
    tx2.error(i) = sqrt((tx2.x_hat(i)-tx2.x).^2 + (tx2.y_hat(i)-tx2.y).^2);
    
    cumulative_class_error(i) = sum(DOA1_hat(1:i)~=DOA1(1:i));
end

%Check if the estimated need swapped
if tx1.error(end) > sqrt((tx1.x_hat(i)-tx2.x).^2 +(tx1.y_hat(i)-tx2.y).^2)
    tmp.doa1 = DOA1_hat;
    DOA1_hat = DOA2_hat;
    DOA2_hat = tmp.doa1;
    
    tmp.y_hat = tx1.y_hat;
    tmp.x_hat = tx1.x_hat;
    
    tx1.y_hat = tx2.y_hat;
    tx1.x_hat = tx2.x_hat;
    tx1.error = sqrt((tx1.x_hat - tx1.x).^2 + (tx1.y_hat - tx1.y).^2);
    
    tx2.y_hat = tmp.y_hat;
    tx2.x_hat = tmp.x_hat; 
    tx2.error = sqrt((tx2.x_hat - tx2.x).^2 + (tx2.y_hat - tx2.y).^2);    
    
    tmp.a = a;
    tmp.b = b;
    
%     a(3,:) = tmp.a(4,:);
%     a(4,:) = tmp.a(3,:);
%     b(3,:) = tmp.b(4,:);
%     b(4,:) = tmp.b(3,:);
    
    tmp.tru = tx1.true_curve;
    
    tx1.true_curve = tx2.true_curve;
    
    tx2.true_curve = tmp.tru;
    
end

figure; 
plot(t, DOA1_hat, 'r');
hold on;
plot(t, DOA2_hat, 'b');
plot(t, fun([a(1,end), b(1,end)], t), 'c-.', 'LineWidth',3);
plot(t, fun([a(2,end), b(2,end)], t), 'c--', 'LineWidth',3);
%plot(t, fun([a(3,end), b(3,end)], t), 'r', 'LineWidth',3);
%plot(t, fun([a(4,end), b(4,end)], t), 'b', 'LineWidth',3);
plot(t, tx1.true_curve);
plot(t, tx2.true_curve);
title('DoA Estimates and Fitted Curves')
xlabel('Time Step')
ylabel('DoA (radians)')
legend('Clustered DOA1 Estimates','Clustered DOA2 Estimates',...
    'Clustered Fitted DOA1','Clustered Fitted DOA2',...
    'True DOA1','True DOA2',...
    'Location','SouthEast');
%'True Fitted DOA1','True Fitted DOA2',...


figure; 
semilogy(t, tx1.error, t, tx2.error);
title('Target Localization Error')
legend('Target One', 'Target Two')

figure
plot(rx.pos_x, rx.pos_y, 'o')
hold on
plot(tx1.x, tx1.y,'ko','MarkerFaceColor','k')
plot(tx2.x, tx2.y,'o','MarkerFaceColor','r')
plot(tx1.x_hat(1:length(t)/10:end), tx1.y_hat(1:length(t)/10:end),'kx')
plot(tx2.x_hat(1:length(t)/10:end), tx2.y_hat(1:length(t)/10:end),'rx')
legend('Collector Position','Target 1 - True Position','Target 2 - True Position', 'Target 1 - Estimated Pos. over Time', 'Target 2 - Estimated Pos. over Time')
title('Position Estimates over Time')

figure; 
plot(t, cumulative_class_error);
title('Cumulative Classification Error');
xlabel('Time');
ylabel('Error');
legend('Cumulative Error');

disp("Final Estimates: ");
fprintf('tx1 x: %f\n', tx1.x_hat(end));
fprintf('tx1 y: %f\n', tx1.y_hat(end));
fprintf('tx2 x: %f\n', tx2.x_hat(end));
fprintf('tx2 y: %f\n', tx2.y_hat(end));
disp("Final Error: ");
fprintf('tx1 error %f\n', tx1.error(end));
fprintf('tx2 error %f\n', tx2.error(end));
disp("Classification Error: ");
fprintf('%2.2f%%\n',(1-sum(DOA1_hat==DOA1)/length(DOA1)));

% %fun = @(x, t)x(1)*atan((t-x(2))/x(3));
% fun = @(x, t)atan((t-x(1))/x(2));
% %initialGuess = [54, 300, 450];
% initialGuess = [100, 100];
% 
% xdata = DOA1;
% x = lsqcurvefit(fun, initialGuess, t(1:end), xdata(1:end));
% ypred = x(1)*Collector.speed + pos.y(1);
% xpred = pos.x(1) - x(2)*Collector.speed;
% pred_err = sqrt(((xpred-Target1.x).^2 +(ypred-Target1.y).^2));
% 
% figure; 
% plot(t, xdata, 'r')
% hold on
% plot(t, fun(x, t), 'b-')

% if saveFlg
%     savename = "data/"+date+"_"+saveStr;
%     save(savename);
% end

if saveFlg
    i=0;
    savename = "data/"+date+"_"+saveStr + "_" + string(i);
    while isfile(savename+".mat")
        i=i+1;
        savename = "data/"+date+"_"+saveStr+"_"+string(i);
    end
    
    save(savename)
    
    display('Data saved as '+savename)
end
    
    
    
    
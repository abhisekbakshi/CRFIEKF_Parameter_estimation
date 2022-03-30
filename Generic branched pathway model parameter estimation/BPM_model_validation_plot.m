function BPM_model_validation_plot(x1, x2)

%%% Normalization
maxx = max(x1);
minn = min(x1);
for i=1:4
    for j=1:5001
       x1_norm(j,i) = (x1(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end

maxx = max(x2);
minn = min(x2);
for i=1:4
    for j=1:5001
       x2_norm(j,i) = (x2(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end



figure
subplot(2,2,1);
plot (x1_norm(:,1),'b')
hold on
plot (x2_norm(:,1),'r')
legend('x1 (CRFIEKF)','x1 (Jia et al. 2012)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(2,2,2);
plot (x1_norm(:,2),'b')
hold on
plot (x2_norm(:,2),'r')
legend('x2 (CRFIEKF)','x2 (Jia et al. 2012)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(2,2,3);
plot (x1_norm(:,3),'b')
hold on
plot (x2_norm(:,3),'r')
legend('x3 (CRFIEKF)','x3 (Jia et al. 2012)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(2,2,4);
plot (x1_norm(:,4),'b')
hold on
plot (x2_norm(:,4),'r')
legend('x4 (CRFIEKF)','x4 (Jia et al. 2012)');
ylabel('Activity/Concentration')
xlabel('Time')


end
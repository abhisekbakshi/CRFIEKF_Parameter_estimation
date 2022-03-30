function IL6_jak_STAT3_validation_plot(x1, x2)

%%% Normalization
maxx = max(x1);
minn = min(x1);
for i=1:6
    for j=1:5001
       x1_norm(j,i) = (x1(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end

maxx = max(x2);
minn = min(x2);
for i=1:6
    for j=1:5001
       x2_norm(j,i) = (x2(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end



figure
subplot(3,2,1);
plot (x1_norm(:,1),'b')
hold on
plot (x2_norm(:,1),'r')
legend('IL-6~IL6R (CRFIEKF)','IL-6~IL6R (Haeseler et al. 2018)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(3,2,2);
plot (x1_norm(:,2),'b')
hold on
plot (x2_norm(:,2),'r')
legend('gp130 (CRFIEKF)','gp130 (Haeseler et al. 2018)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(3,2,3);
plot (x1_norm(:,3),'b')
hold on
plot (x2_norm(:,3),'r')
legend('(p)Rcomplex (CRFIEKF)','(p)Rcomplex (Haeseler et al. 2018)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(3,2,4);
plot (x1_norm(:,4),'b')
hold on
plot (x2_norm(:,4),'r')
legend('(p)STAT3 (CRFIEKF)','(p)STAT3 (Haeseler et al. 2018)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(3,2,5);
plot (x1_norm(:,5),'b')
hold on
plot (x2_norm(:,5),'r')
legend('SOCS3 mRNA (CRFIEKF)','SOCS3 mRNA (Haeseler et al. 2018)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(3,2,6);
plot (x1_norm(:,6),'b')
hold on
plot (x2_norm(:,6),'r')
legend('SOCS3 protein (CRFIEKF)','SOCS3 protein (Haeseler et al. 2018)');
ylabel('Activity/Concentration')
xlabel('Time')


end
function Parameter_estimation_FIS_membership_vs_Kinoshita_plot(x1, x2, x3, x4)

%%% Extracts dynamics of 8 metabolites 
x1 = [x1(:,1) x1(:,3) x1(:,4) x1(:,6) x1(:,9) x1(:,11) x1(:,12) x1(:,13)];
x2 = [x2(:,1) x2(:,3) x2(:,4) x2(:,6) x2(:,9) x2(:,11) x2(:,12) x2(:,13)];
x3 = [x3(:,1) x3(:,3) x3(:,4) x3(:,6) x3(:,9) x3(:,11) x3(:,12) x3(:,13)];
x4 = [x4(:,1) x4(:,3) x4(:,4) x4(:,6) x4(:,9) x4(:,11) x4(:,12) x4(:,13)];



%%% Normalization
maxx = max(x1);
minn = min(x1);
for i=1:8
    for j=1:5001
       x1_norm(j,i) = (x1(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end

maxx = max(x2);
minn = min(x2);
for i=1:8
    for j=1:5001
       x2_norm(j,i) = (x2(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end

maxx = max(x3);
minn = min(x3);
for i=1:8
    for j=1:5001
       x3_norm(j,i) = (x3(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end

maxx = max(x4);
minn = min(x4);
for i=1:8
    for j=1:5001
       x4_norm(j,i) = (x4(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end


%%% Load Kinoshita simulation Result
Kinoshita_validation = load('Kinoshita_simulation_Result.mat');
for v = fieldnames(Kinoshita_validation)
    hypoxia_norm = Kinoshita_validation.(v{1});
end





figure
subplot(4,2,1);
plot (x1_norm(:,1),'b')
hold on
plot (x2_norm(:,1),'k')
hold on
plot (x3_norm(:,1),'r')
hold on
plot (x4_norm(:,1),'m')
hold on
plot (hypoxia_norm(:,1),'g')
legend('Glucose 6P with Gaussian','Glucose 6P with Generalised Bell','Glucose 6P with Triangular','Glucose 6P with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,2);
plot (x1_norm(:,2),'b')
hold on
plot (x2_norm(:,2),'k')
hold on
plot (x3_norm(:,2),'r')
hold on
plot (x4_norm(:,2),'m')
hold on
plot (hypoxia_norm(:,2),'g')
legend('Fructose 6P with Gaussian','Fructose 6P with Generalised Bell','Fructose 6P with Triangular','Fructose 6P with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,3);
plot (x1_norm(:,3),'b')
hold on
plot (x2_norm(:,3),'k')
hold on
plot (x3_norm(:,3),'r')
hold on
plot (x4_norm(:,3),'m')
hold on
plot (hypoxia_norm(:,3),'g')
legend('Fructose 1,6BP with Gaussian','Fructose 1,6BP with Generalised Bell','Fructose 1,6BP with Triangular','Fructose 1,6BP with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,4);
plot (x1_norm(:,4),'b')
hold on
plot (x2_norm(:,4),'k')
hold on
plot (x3_norm(:,4),'r')
hold on
plot (x4_norm(:,4),'m')
hold on
plot (hypoxia_norm(:,4),'g')
legend('DHAP with Gaussian','DHAP with Generalised Bell','DHAP with Triangular','DHAP with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x1_norm(:,5),'b')
hold on
plot (x2_norm(:,5),'k')
hold on
plot (x3_norm(:,5),'r')
hold on
plot (x4_norm(:,5),'m')
hold on

plot (hypoxia_norm(:,5),'g')
legend('3PG with Gaussian','3PG with Generalised Bell','3PG with Triangular','3PG with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,6);
plot (x1_norm(:,6),'b')
hold on
plot (x2_norm(:,6),'k')
hold on
plot (x3_norm(:,6),'r')
hold on
plot (x4_norm(:,6),'m')
hold on
plot (hypoxia_norm(:,6),'g')
legend('PEP with Gaussian','PEP with Generalised Bell','PEP with Triangular','PEP with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x1_norm(:,7),'b')
hold on
plot (x2_norm(:,7),'k')
hold on
plot (x3_norm(:,7),'r')
hold on
plot (x4_norm(:,7),'m')
hold on
plot (hypoxia_norm(:,7),'g')
legend('Pyruvate with Gaussian','Pyruvate with Generalised Bell','Pyruvate with Triangular','Pyruvate with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,8);
plot (x1_norm(:,8),'b')
hold on
plot (x2_norm(:,8),'k')
hold on
plot (x3_norm(:,8),'r')
hold on
plot (x4_norm(:,8),'m')
hold on
plot (hypoxia_norm(:,8),'g')
legend('Lactate with Gaussian','Lactate with Generalised Bell','Lactate with Triangular','Lactate with Trapiziodal','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')

end


function Parameter_estimation_by_6_vs_Kinoshita_vs_plot(x1, x2)

%%% Extracts dynamics of 8 metabolites 
x1 = [x1(:,1) x1(:,3) x1(:,4) x1(:,6) x1(:,9) x1(:,11) x1(:,12) x1(:,13)];
x2 = [x2(:,1) x2(:,3) x2(:,4) x2(:,6) x2(:,9) x2(:,11) x2(:,12) x2(:,13)];


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


%%% Load Kinoshita simulation Result
Kinoshita_validation = load('Kinoshita_simulation_Result.mat');
for v = fieldnames(Kinoshita_validation)
    hypoxia_norm = Kinoshita_validation.(v{1});
end


%%% Load CE MS measurement original simulation Result
CE_MS_dataset = load('CE_MS_Measurement_Original_dataset.mat');
for v = fieldnames(CE_MS_dataset)
    CE_MS_data = CE_MS_dataset.(v{1});
end







figure
subplot(4,2,1);
plot (x1_norm(:,1),'b')
hold on
plot (x2_norm(:,1),'r')
hold on
plot (hypoxia_norm(:,1),'g')
legend('Glucose 6P with CRFIEKF','Glucose 6P with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,2);
plot (CE_MS_data(:,1),'k')
legend('Glucose 6P CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,3);
plot (x1_norm(:,2),'b')
hold on
plot (x2_norm(:,2),'r')
hold on
plot (hypoxia_norm(:,2),'g')
legend('Fructose 6P with CRFIEKF','Fructose 6P with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,4);
plot (CE_MS_data(:,2),'k')
legend('Fructose 6P CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x1_norm(:,3),'b')
hold on
plot (x2_norm(:,3),'r')
hold on
plot (hypoxia_norm(:,3),'g')
legend('Fructose 1,6BP with CRFIEKF','Fructose 1,6BP with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,6);
plot (CE_MS_data(:,3),'k')
legend('Fructose 1,6BP CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x1_norm(:,4),'b')
hold on
plot (x2_norm(:,4),'r')
hold on
plot (hypoxia_norm(:,4),'g')
legend('DHAP with CRFIEKF','DHAP with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,8);
plot (CE_MS_data(:,4),'k')
legend('DHAP CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')




figure
subplot(4,2,1);
plot (x1_norm(:,5),'b')
hold on
plot (x2_norm(:,5),'r')
hold on
plot (hypoxia_norm(:,5),'g')
legend('3PG with CRFIEKF','3PG with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,2);
plot (CE_MS_data(:,5),'k')
legend('3PG CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,3);
plot (x1_norm(:,6),'b')
hold on
plot (x2_norm(:,6),'r')
hold on
plot (hypoxia_norm(:,6),'g')
legend('PEP with CRFIEKF','PEP with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,4);
plot (CE_MS_data(:,6),'k')
legend('PEP CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x1_norm(:,7),'b')
hold on
plot (x2_norm(:,7),'r')
hold on
plot (hypoxia_norm(:,7),'g')
legend('Pyruvate with CRFIEKF','Pyruvate with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,6);
plot (CE_MS_data(:,7),'k')
legend('Pyruvate CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x1_norm(:,8),'b')
hold on
plot (x2_norm(:,8),'r')
hold on
plot (hypoxia_norm(:,8),'g')
legend('Lactate with CRFIEKF','Lactate with HEKF','Kinoshita et. al');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,8);
plot (CE_MS_data(:,8),'k')
legend('Lactate CE MS (Kinoshita et. al)');
ylabel('Activity/Concentration')
xlabel('Time')






% % % 
% % % subplot(4,2,3);
% % % plot (x1_norm(:,3),'b')
% % % hold on
% % % plot (x2_norm(:,3),'r')
% % % hold on
% % % plot (hypoxia_norm(:,3),'g')
% % % legend('Fructose 1,6BP with 6 FIS','Fructose 1,6BP with 5 FIS','Fructose 1,6BP with 4 FIS','Kinoshita et. al');
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % subplot(4,2,4);
% % % plot (x1_norm(:,4),'b')
% % % hold on
% % % plot (x2_norm(:,4),'k')
% % % hold on
% % % plot (x3_norm(:,4),'r')
% % % hold on
% % % plot (hypoxia_norm(:,4),'g')
% % % legend('DHAP with 6 FIS','DHAP with 5 FIS','DHAP with 4 FIS','Kinoshita et. al');
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % subplot(4,2,5);
% % % plot (x1_norm(:,5),'b')
% % % hold on
% % % plot (x2_norm(:,5),'k')
% % % hold on
% % % plot (x3_norm(:,5),'r')
% % % hold on
% % % plot (hypoxia_norm(:,5),'g')
% % % legend('3PG with 6 FIS','3PG with 5 FIS','3PG with 4 FIS','Kinoshita et. al');
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % subplot(4,2,6);
% % % plot (x1_norm(:,6),'b')
% % % hold on
% % % plot (x2_norm(:,6),'k')
% % % hold on
% % % plot (x3_norm(:,6),'r')
% % % hold on
% % % plot (hypoxia_norm(:,6),'g')
% % % legend('PEP with 6 FIS','PEP with 5 FIS','PEP with 4 FIS','Kinoshita et. al');
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % subplot(4,2,7);
% % % plot (x1_norm(:,7),'b')
% % % hold on
% % % plot (x2_norm(:,7),'k')
% % % hold on
% % % plot (x3_norm(:,7),'r')
% % % hold on
% % % plot (hypoxia_norm(:,7),'g')
% % % legend('Pyruvate with 6 FIS','Pyruvate with 5 FIS','Pyruvate with 4 FIS','Kinoshita et. al');
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')
% % % 
% % % 
% % % subplot(4,2,8);
% % % plot (x1_norm(:,8),'b')
% % % hold on
% % % plot (x2_norm(:,8),'k')
% % % hold on
% % % plot (x3_norm(:,8),'r')
% % % hold on
% % % plot (hypoxia_norm(:,8),'g')
% % % legend('Lactate with 6 FIS','Lactate with 5 FIS','Lactate with 4 FIS','Kinoshita et. al');
% % % ylabel('Activity/Concentration')
% % % xlabel('Time')

end


function [] = JA_STAT_model_validation_plot(x1, total_stat5, ph_stat5)


%%% Normalization
maxx = max(x1);
minn = min(x1);
for i=1:5
    for j=1:16
       x_v(j,i) = (x1(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end



%%% Load Waightlab dataset from Timmer et al. 20004
Weightlab_original = load('Weightlab_original_dataset.mat');
for v = fieldnames(Weightlab_original)
    original = Weightlab_original.(v{1});
end

%%% Load molecular dynamics presented in Sun et al. (2008)
Molecule_activity_original = load('Molecule_activity_plosone_model.mat');
for v = fieldnames(Molecule_activity_original)
    plosOne = Molecule_activity_original.(v{1});
end

%%% Load molecular dynamics presented in Sun et al. (2008) and Sun et al. (2016)
Molecule_activity_previous_model = load('Molecule_activity_previous_model.mat');
for v = fieldnames(Molecule_activity_previous_model)
    Other_model_normalized = Molecule_activity_previous_model.(v{1});
    Sun_2016_normalized = Molecule_activity_previous_model.(v{2});
end




figure
subplot(2,2,1);
plot (x_v(:,1),'b')
hold on
plot (Other_model_normalized(:,1),'r')
hold on
plot (Sun_2016_normalized(:,1),'k')
legend('STAT5 (CRFIEKF)','STAT5 (Sun et. al 2008)','STAT5 (Sun et. al 2016)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(2,2,2);
plot (x_v(:,2),'b')
hold on
plot (Other_model_normalized(:,2),'r')
hold on
plot (Sun_2016_normalized(:,2),'k')
legend('STAT5P (CRFIEKF)','STAT5P (Sun et. al 2008)','STAT5P (Sun et. al 2016)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(2,2,3);
plot (x_v(:,3),'b')
hold on
plot (Other_model_normalized(:,3),'r')
hold on
plot (Sun_2016_normalized(:,3),'k')
legend('STAT5 2P (CRFIEKF)','STAT5 2P (Sun et. al 2008)','STAT5 2SP (Sun et. al 2016)');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(2,2,4);
plot (x_v(:,4),'b')
hold on
plot (Other_model_normalized(:,4),'r')
hold on
plot (Sun_2016_normalized(:,4),'k')
legend('STAT5N (CRFIEKF)','STAT5N (Sun et. al 2008)','STAT5N (Sun et. al 2016)');
ylabel('Activity/Concentration')
xlabel('Time')




figure
subplot(2,1,1);
plot (total_stat5(:,1),'b')
hold on
plot (original(:,2),'k')
hold on
plot (plosOne(:,2),'r')
legend('Total STAT5 (CRFIEKF)','Total STAT5 (Original)','Total STAT5 (PlosOne)');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(2,1,2);
plot (ph_stat5(:,1),'b')
hold on
plot (original(:,1),'k')
hold on
plot (plosOne(:,1),'r')
legend('Tyrosine STAT5 (CRFIEKF)','Tyrosine STAT5 (Original)','Tyrosine STAT5 (PlosOne)');
ylabel('Activity/Concentration')
xlabel('Time')


end
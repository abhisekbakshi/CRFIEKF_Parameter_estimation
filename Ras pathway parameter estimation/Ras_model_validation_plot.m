function [x_v] = Ras_model_validation_plot(x_v)

%%% Load weightlab data Sun et al. (2008)
Ras_validation = load('Ras_weight_lab_data.mat');
for v = fieldnames(Ras_validation)
    Ras_x_v = Ras_validation.(v{1});
end


figure
subplot(2,1,1);
plot (Ras_x_v(:,1),'b')
hold on;
plot (x_v(:,1),'r')
legend('Raf-1* observed concentration','Raf-1* predicted concentration');
ylabel('Activity/Concentration')
xlabel('Time')

subplot(2,1,2);
plot (Ras_x_v(:,2),'b')
hold on;
plot (x_v(:,2),'r')
legend('RKIP observed concentration','RKIP predicted concentration');
ylabel('Activity/Concentration')
xlabel('Time')

end



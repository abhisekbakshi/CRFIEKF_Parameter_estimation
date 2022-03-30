function [] = Yeast_glyolysis_validation_plot(x_v, reference, reference_plosone_noiseless, reference_plosone_noisy)


%%% Normalization
maxx = max(x_v);
minn = min(x_v);
for i=1:7
    for j=1:5001
       x_v_norm(j,i) = (x_v(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end


maxx = max(reference);
minn = min(reference);
for i=1:7
    for j=1:5001
       reference_norm(j,i) = (reference(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end


maxx = max(reference_plosone_noiseless);
minn = min(reference_plosone_noiseless);
for i=1:7
    for j=1:5001
       reference_plosone_noiseless_norm(j,i) = (reference_plosone_noiseless(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end


maxx = max(reference_plosone_noisy);
minn = min(reference_plosone_noisy);
for i=1:7
    for j=1:5001
       reference_plosone_noisy_norm(j,i) = (reference_plosone_noisy(j,i) - minn(1,i))/(maxx(1,i) - minn(1,i));
    end
end


figure
subplot(4,2,1);
plot (x_v_norm(:,1),'b')
hold on
plot (reference_norm(:,1),'m')
hold on
plot (reference_plosone_noiseless_norm(:,1),'k')
hold on
plot (reference_plosone_noisy_norm(:,1),'r')
legend('Glucose our model','Glucose (Ruoff et al. (2003))','Glucose (Yazdani et al. (2020))','Glucose (Yazdani et al. (2020))');
ylabel('Activity/Concentration')
xlabel('Time')




subplot(4,2,2);
plot (x_v_norm(:,2),'b')
hold on
plot (reference_norm(:,2),'m')
hold on
plot (reference_plosone_noiseless_norm(:,2),'k')
hold on
plot (reference_plosone_noisy_norm(:,2),'r')
legend('Glyceraldehyde 3P our model','Glyceraldehyde 3P (Ruoff et al. (2003))','Glyceraldehyde 3P (Yazdani et al. (2020))','Glyceraldehyde 3P (Yazdani et al. (2020))');
ylabel('Activity/Concentration')
xlabel('Time')



subplot(4,2,3);
plot (x_v_norm(:,3),'b')
hold on
plot (reference_norm(:,3),'m')
hold on
plot (reference_plosone_noiseless_norm(:,3),'k')
hold on
plot(reference_plosone_noisy_norm(:,3),'r')
legend('BPG 1,3 our model','BPG 1,3 (Ruoff et al. (2003))','BPG 1,3 (Yazdani et al. (2020))','BPG 1,3 (Yazdani et al. (2020))');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,4);
plot (x_v_norm(:,4),'b')
hold on
plot (reference_norm(:,4),'m')
hold on
plot (reference_plosone_noiseless_norm(:,4),'k')
hold on
plot (reference_plosone_noisy_norm(:,4),'r')
legend('Pyruvate our model','Pyruvate (Ruoff et al. (2003))','Pyruvate (Yazdani et al. (2020))','Pyruvate (Yazdani et al. (2020))');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,5);
plot (x_v_norm(:,5),'b')
hold on
plot (reference_norm(:,5),'m')
hold on
plot (reference_plosone_noiseless_norm(:,5),'k')
hold on
plot (reference_plosone_noisy_norm(:,5),'r')
legend('NADH our model','NADH (Ruoff et al. (2003))','NADH (Yazdani et al. (2020))','NADH (Yazdani et al. (2020))');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,6);
plot (x_v_norm(:,6),'b')
hold on
plot (reference_norm(:,6),'m')
hold on
plot (reference_plosone_noiseless_norm(:,6),'k')
hold on
plot (reference_plosone_noisy_norm(:,6),'r')
legend('ATP our model','ATP (Ruoff et al. (2003))','ATP (Yazdani et al. (2020))','ATP (Yazdani et al. (2020))');
ylabel('Activity/Concentration')
xlabel('Time')


subplot(4,2,7);
plot (x_v_norm(:,7),'b')
hold on
plot (reference_norm(:,7),'m')
hold on
plot (reference_plosone_noiseless_norm(:,7),'k')
hold on
plot (reference_plosone_noisy_norm(:,7),'r')
legend('Pyruvate ext. our model','Pyruvate ext. (Ruoff et al. (2003))','Pyruvate ext. (Yazdani et al. (2020))','Pyruvate ext. (Yazdani et al. (2020))');
ylabel('Activity/Concentration')
xlabel('Time')



end


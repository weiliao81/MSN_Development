clc;clear;
load 'gradient1-3 data'
load 'subject information'
nsubs = 1790;
network = importdata('von Economo index');
withinDispersion = zeros(1790,1533);
withinDispersion_mean = zeros(7,1790);
%%
%**********************************************step1 : Caculate Within Dispersion
    for i = 1:nsubs
        for j=1:7
            index = find(network == j);
            cen1 = median(gradientAll1(i,index));
            cen2 = median(gradientAll2(i,index));
            cen3 = median(gradientAll3(i,index));
            for k=1:length(index)
                distance(index(k),i) = sqrt((gradientAll1(i,index(k))- cen1).^2 + ...
                    (gradientAll2(i,index(k))- cen2).^2 + ...
                    (gradientAll3(i,index(k))- cen3).^2 );
            end
            withinDispersion_mean(j,i) = mean(distance(index,i));
        end
    end
    mod = [sex eTIV];
    withinDispersion_harmonized = combat(withinDispersion_mean, batch', mod,1);
%%
%**********************************************step2 : Save Data
withinResult = withDispersion_harmonized'
withinNames = {'age','eTIV','sex','	motor','asso1','asso2','secSensory',...
    'primSensory','limbic','insula'};
data = array2table([age sex eTIV withinResult],'VariableNames',withinNames);
writetable(data, 'withinDispersion.csv')


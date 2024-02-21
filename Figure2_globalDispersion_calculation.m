clc;clear;
load 'gradient1-3 data'
load 'subject infomation'
nsubs = 1790;
nROIs = 1533;
%%
%**********************************************step1 : Caculate Global Dispersion

    for i=1:nsubs
        
        cen1 = median(gradientAll1(i,:));
        cen2 = median(gradientAll2(i,:));
        cen3 = median(gradientAll3(i,:));
        
        for j=1:nROIs
            distance(j,i) = sqrt((gradientAll1(i,j)- cen1).^2 + ...
                (gradientAll2(i,j)- cen2).^2 + ...
                (gradientAll3(i,j)- cen3).^2 );
        end
        
    end
    mod = [sex eTIV];
    distance_harmonized = combat(distance, batch', mod,1);

%%
%**********************************************step2 : Save Data
columes_name = {'age','eTIV','sex','globalDispersion'};
result_global = table(age ,eTIV ,sex ,mean(distance_harmonized)',...
    'VariableNames',columes_name);
writetable(result_global,'global dispersion.csv')

result_wholebrain = table(age ,eTIV ,sex ,(distance_harmonized)',...
    'VariableNames',columes_name);
writetable(result_wholebrain,'whole brain regions' dispersion.csv')
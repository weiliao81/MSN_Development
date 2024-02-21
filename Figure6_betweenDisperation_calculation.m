
clc;clear;
load 'gradient1-3 data'
load 'subjects information'
mod = [sex eTIV];
von_parc = load('vonEconomo index');
%%
%**********************************************step1 : Caculate Between Dispersion
    for i=1:1790
        for k=1:7
            index = find(von_parc == k);
            gradientCen1(i,k) = median(gradientAll1(i,index));
            gradientCen2(i,k) = median(gradientAll2(i,index));
            gradientCen3(i,k) = median(gradientAll3(i,index));
        end
        clear index;
    end
    
    von_cog(:,1,:) = gradientCen1;
    von_cog(:,2,:) = gradientCen2;
    von_cog(:,3,:) = gradientCen3;
    net_dist = zeros(7,7,length(age));
  
    for n1 = 1:7
        for n2 = 1:7
            net_dist(n1,n2,:) = sqrt(sum((von_cog(:,:,n1)' - von_cog(:,:,n2)').^2));
        end
    end
    for s = 1:length(age)
        net_dist_long(s,:) = squareform(net_dist(:,:,s));
    end
    net_dist_long_harm = combat(net_dist_long', batch', mod,1);
    
%%
%**********************************************step2 : Save Data
beResult = net_dist_long_harm';
between_names ={'age','sex','eTIV',...
    'b12','b13','b14','b15','b16','b17','b23'...
    'b24','b25','b26','b27','b34','b35','b36'...
    'b37','b45','b46','b47','b56','b57','b67'};
data = array2table([age sex eTIV beResult],'VariableNames',between_names);
writetable(data, 'betweenDispersion.csv')


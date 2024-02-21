clc;clear all;close all;
%%%% this code copied from https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md

Tvalue=csvread('whole brain regions' gamlssTvalues.csv');
gene_data=csvread('brain regions * gene expression.csv',1,1);
y=Tvalue(1:1:767);
x=gene_data.data;
x(768:1:end,:)=[];
x=zscore(x);
y=zscore(y);
genes=gene_data.textdata;
geneindex=1:size(x,2);
%%
%**********************************************step1 : caculate pls scores
for f2_plsScores=1
    dim=2;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(x,y,dim);
    [R1,p1]=corr([XS(:,1),XS(:,2)],y);
    if R1(1,1)<0
        stats.W(:,1)=-1*stats.W(:,1);
        XS(:,1)=-1*XS(:,1);
    end
    if R1(2,1)<0
        stats.W(:,2)=-1*stats.W(:,2);
        XS(:,2)=-1*XS(:,2);
    end
    %print out results 
    csvwrite('data/gene/PLS1_ROIscores.csv',XS(:,1));
    csvwrite('data/gene/PLS2_ROIscores.csv',XS(:,2));
end

%%
%**********************************************step2 : caculate pls weights
for f2_plsWeights=1
    [PLS1w,x1] = sort(stats.W(:,1),'descend');
    PLS1ids=genes(x1);
    geneindex1=geneindex(x1);
    [PLS2w,x2] = sort(stats.W(:,2),'descend');
    PLS2ids=genes(x2);
    geneindex2=geneindex(x2);
    
    %define variables for storing the (ordered) weights from all bootstrap runs
    PLS1weights=[];
    PLS2weights=[];
    
    %start bootstrap
    bootnum=1000
    for i=1:bootnum
        i
        myresample = randsample(size(x,1),size(x,1));
        res(i,:)=myresample; %store resampling out of interest
        Xr=x(myresample,:); % define X for resampled subjects
        Yr=y(myresample,:); % define X for resampled subjects
        dim=2;
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data 
        temp=stats.W(:,1);%extract PLS1 weights
        newW=temp(x1); %order the newly obtained weights the same way as initial PLS
        if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW=-1*newW;
        end
        PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
        temp=stats.W(:,2);%extract PLS2 weights
        newW=temp(x2); %order the newly obtained weights the same way as initial PLS
        if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW=-1*newW;
        end
        PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run
    end
    %get standard deviation of weights from bootstrap runs
    PLS1sw=std(PLS1weights');
    PLS2sw=std(PLS2weights');
    %get bootstrap weights
    temp1=PLS1w./PLS1sw';
    temp2=PLS2w./PLS2sw';
    %order bootstrap weights (Z) and names of regions
    [Z1 ind1]=sort(temp1,'descend');
    PLS1=PLS1ids(ind1);
    geneindex1=geneindex1(ind1);
    [Z2 ind2]=sort(temp2,'descend');
    PLS2=PLS2ids(ind2);
    geneindex2=geneindex2(ind2);
    
    %print out results
    fid1 = fopen(['data/gene/PLS1_geneWeights_','1000','.csv'],'w')
    for i=1:length(genes)
        fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
    end
    fclose(fid1)
    
    fid2 = fopen(['data/gene/PLS2_geneWeights_','1000','.csv'],'w')
    for i=1:length(genes)
        fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
    end
    fclose(fid2)
end

%COMPARETRACES Function outputs parameters that compare the automated trace
%to the gold standard trace.

function IndexStr = CompareTraces(AM_G, r_G, R_G, AM_A, r_A, R_A, h, delta)

% [AM_G r_G R_G ~] = Eliminate_Terminal_Branches(AM_G,r_G,R_G,h,1,0); 
% [AM_A r_A R_A ~] = Eliminate_Terminal_Branches(AM_A,r_A,R_A,h,1,0);

[AM_G r_G R_G] = AdjustPPM(AM_G,r_G,R_G,1/delta);
[AM_A r_A R_A] = AdjustPPM(AM_A,r_A,R_A,1/delta);

temp=Compare2Trees(AM_G,r_G,R_G,AM_A,r_A,R_A,h,delta);

AMlbl_trees_G = LabelTreesAM(AM_G);
Labels_G=unique(AMlbl_trees_G(AMlbl_trees_G>0));

AMlbl_trees_A = LabelTreesAM(AM_A);
Labels_A=unique(AMlbl_trees_A(AMlbl_trees_A>0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters of individual traces 
IndexStr.N_trees_G=length(Labels_G);
IndexStr.N_trees_A=length(Labels_A);
IndexStr.Lg=temp.Lg;
IndexStr.La=temp.La;
IndexStr.BPg=temp.BPg;
IndexStr.BPa=temp.BPa;
IndexStr.TPg=temp.TPg;
IndexStr.TPa=temp.TPa;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of entire traces
IndexStr.corresponding_length_ing=temp.corresponding_length_ing;
IndexStr.corresponding_length_ina=temp.corresponding_length_ina;
IndexStr.Tortuosity_ag=IndexStr.corresponding_length_ina/IndexStr.corresponding_length_ing-1;
IndexStr.Dag_trace_mean=temp.Dag_trace_mean;
IndexStr.Dag_trace_SD=temp.Dag_trace_SD;
IndexStr.fp_length=temp.fp_length;
IndexStr.fn_length=temp.fn_length;
IndexStr.corresponding_BP_ing=temp.corresponding_BP_ing;
IndexStr.corresponding_BP_ina=temp.corresponding_BP_ina;
IndexStr.Dag_BP_mean=temp.Dag_BP_mean;
IndexStr.Dag_BP_SD=temp.Dag_BP_SD;
IndexStr.fp_BP=temp.fp_BP;
IndexStr.fn_BP=temp.fn_BP;
IndexStr.corresponding_TP_ing=temp.corresponding_TP_ing;
IndexStr.corresponding_TP_ina=temp.corresponding_TP_ina;
IndexStr.Dag_TP_mean=temp.Dag_TP_mean;
IndexStr.Dag_TP_SD=temp.Dag_TP_SD;
IndexStr.fp_TP=temp.fp_TP;
IndexStr.fn_TP=temp.fn_TP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of corresponding trees

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of corresponding branches

% % NaN_temp=nan(IndexStr.N_trees_G,IndexStr.N_trees_A);
% % TracesIndexStr=struct('Lg',NaN_temp,'La',NaN_temp,'BPg',NaN_temp,'BPa',NaN_temp, ...
% %     'TPg',NaN_temp,'TPa',NaN_temp,'Dag_trace_mean',NaN_temp,'Dag_trace_SD',NaN_temp, ...
% %     'corresponding_length_ing',NaN_temp,'corresponding_length_ina',NaN_temp, ...
% %     'fp_length',NaN_temp,'fn_length',NaN_temp,'fn_TP',NaN_temp,'fp_TP',NaN_temp, ...
% %     'Dag_TP_mean',NaN_temp,'Dag_TP_SD',NaN_temp,'fn_BP',NaN_temp,'fp_BP',NaN_temp, ...
% %     'Dag_BP_mean',NaN_temp,'Dag_BP_SD',NaN_temp);
% % for i=1:length(Labels_G)
% %     [e1,e2]=find(AMlbl_trees_G==Labels_G(i));
% %     ee=unique([e1;e2]);
% %     AMg=AMlbl_trees_G(ee,ee);
% %     %AMg(AMg>0)=1; 
% %     rg=r_G(ee,:);
% %     Rg=R_G(ee);
% %     for j=1:length(Labels_A), [i,j] 
% %         [e1,e2]=find(AMlbl_trees_A==Labels_A(j));
% %         ee=unique([e1;e2]);
% %         AMa=AMlbl_trees_A(ee,ee);
% %         %AMa(AMa>0)=1;
% %         ra=r_A(ee,:);
% %         Ra=R_A(ee);
% %         
% %         temp=Compare2Trees(AMg, rg, Rg, AMa, ra, Ra, h, delta);
% %         TracesIndexStr.Lg(i,j)=temp.Lg;
% %         TracesIndexStr.La(i,j)=temp.La;
% %         TracesIndexStr.BPg(i,j)=temp.BPg;
% %         TracesIndexStr.BPa(i,j)=temp.BPa;
% %         TracesIndexStr.TPg(i,j)=temp.TPg;
% %         TracesIndexStr.TPa(i,j)=temp.TPa;
% %         TracesIndexStr.Dag_trace_mean(i,j)=temp.Dag_trace_mean;
% %         TracesIndexStr.Dag_trace_SD(i,j)=temp.Dag_trace_SD;
% %         TracesIndexStr.corresponding_length_ing(i,j)=temp.corresponding_length_ing;
% %         TracesIndexStr.corresponding_length_ina(i,j)=temp.corresponding_length_ina;
% %         TracesIndexStr.fp_length(i,j)=temp.fp_length;
% %         TracesIndexStr.fn_length(i,j)=temp.fn_length;
% %         TracesIndexStr.fn_TP(i,j)=temp.fn_TP;
% %         TracesIndexStr.fp_TP(i,j)=temp.fp_TP;
% %         TracesIndexStr.Dag_TP_mean(i,j)=temp.Dag_TP_mean;
% %         TracesIndexStr.Dag_TP_SD(i,j)=temp.Dag_TP_SD;
% %         TracesIndexStr.fn_BP(i,j)=temp.fn_BP;
% %         TracesIndexStr.fp_BP(i,j)=temp.fp_BP;
% %         TracesIndexStr.Dag_BP_mean(i,j)=temp.Dag_BP_mean;
% %         TracesIndexStr.Dag_BP_SD(i,j)=temp.Dag_BP_SD;
% %     end
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % tree to tree comparisons
% % IndexStr.fp_length=0;
% % IndexStr.fn_length=0;
% % IndexStr.fn_TP=0;
% % IndexStr.fp_TP=0;
% % 
% % IndexStr.fn_BP=0;
% % IndexStr.fp_BP=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IndexStr = Compare2Trees(AMg, rg, Rg, AMa, ra, Ra, h, delta)

%removing AM labels
AMg(AMg>0)=1; 
AMa(AMa>0)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating total length of the circuit
[ig,jg]=find(AMg); 
lllg=AMg; 
lllg(AMg>0)=sum((rg(ig,:)-rg(jg,:)).^2,2).^0.5;
lg=sum(lllg,1)./2;
IndexStr.Lg=full(sum(lg));

[ia,ja]=find(AMa);
llla=AMa; llla(AMa>0)=sum((ra(ia,:)-ra(ja,:)).^2,2).^0.5;
la=sum(llla,2)./2;
IndexStr.La=full(sum(la));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding Branch Points
BPg_ind=(sum(AMg)>2);
IndexStr.BPg=full(sum(BPg_ind));

BPa_ind=(sum(AMa)>2);
IndexStr.BPa=full(sum(BPa_ind));

%Finding Terminal Points
TPg_ind=(sum(AMg)==1);
IndexStr.TPg=full(sum(TPg_ind));

TPa_ind=(sum(AMa)==1);
IndexStr.TPa=full(sum(TPa_ind));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distance matrix calculation
Dag_min=sparse(size(AMg,1),size(AMa,1)); 
temp_length=(size(AMg,1)+size(AMa,1))*ceil(2*h/delta);
Dag_i=zeros(temp_length,1);
Dag_j=zeros(temp_length,1);
Dag_v=zeros(temp_length,1);
count=1;

for j=1:size(AMg,1) 
    d=((ra(:,1)-rg(j,1)).^2+(ra(:,2)-rg(j,2)).^2+(ra(:,3)-rg(j,3)).^2).^0.5;
    ind=find(d<=h);
    if ~isempty(ind)
        Dag_i(count:count+length(ind)-1)=j.*ones(length(ind),1);
        Dag_j(count:count+length(ind)-1)=ind;
        Dag_v(count:count+length(ind)-1)=d(ind);
        count=count+length(ind);
        
        % Calculation of distance from a to g
        [d_min,ind_min]=min(d); 
        Dag_min(j,ind_min)=d_min;     
    end
end
Dag=sparse(Dag_i(1:count-1),Dag_j(1:count-1),Dag_v(1:count-1),size(AMg,1),size(AMa,1)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of distance and false positive/negative lengths
IndexStr.Dag_trace_mean=full((lg*Dag_min*la)/(lg*(Dag_min>0)*la));
IndexStr.Dag_trace_SD=full(((lg*Dag_min.^2*la)/(lg*(Dag_min>0)*la)-IndexStr.Dag_trace_mean^2)^0.5);
IndexStr.corresponding_length_ing=full(lg*(sum(Dag,2)>0));
IndexStr.corresponding_length_ina=full((sum(Dag,1)>0)*la);
IndexStr.fp_length=full((sum(Dag,1)==0)*la);
IndexStr.fn_length=full(lg*(sum(Dag,2)==0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Calculation of distance and numbers of false positive/negative branch points 
Dag_BP=Dag(BPg_ind,BPa_ind);
if ~isempty(Dag_BP)
    IndexStr.corresponding_BP_ing=full(sum(sum(Dag_BP,2)>0));
    IndexStr.corresponding_BP_ina=full(sum(sum(Dag_BP,1)>0));
    IndexStr.fp_BP=full(sum(sum(Dag_BP,1)==0));
    IndexStr.fn_BP=full(sum(sum(Dag_BP,2)==0));
else
    IndexStr.corresponding_BP_ing=0;
    IndexStr.corresponding_BP_ina=0;
    IndexStr.fp_BP=full(sum(BPa_ind));
    IndexStr.fn_BP=full(sum(BPg_ind));
end

Dag_BP_temp=Dag_BP;
Dag_BP_temp(Dag_BP_temp==0)=inf;
Dag_BP_temp=min(Dag_BP_temp,[],2);
Dag_BP_temp(Dag_BP_temp==inf)=[];
if ~isempty(Dag_BP_temp)
    IndexStr.Dag_BP_mean=full(mean(Dag_BP_temp));
    IndexStr.Dag_BP_SD=full(std(Dag_BP_temp));
else
    IndexStr.Dag_BP_mean=NaN;
    IndexStr.Dag_BP_SD=NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of distance and numbers of false positive/negative terminal points 
Dag_TP=Dag(TPg_ind,TPa_ind);
if ~isempty(Dag_TP)
    IndexStr.corresponding_TP_ing=full(sum(sum(Dag_TP,2)>0));
    IndexStr.corresponding_TP_ina=full(sum(sum(Dag_TP,1)>0));
    IndexStr.fp_TP=full(sum(sum(Dag_TP,1)==0));
    IndexStr.fn_TP=full(sum(sum(Dag_TP,2)==0));
else
    IndexStr.corresponding_TP_ing=0;
    IndexStr.corresponding_TP_ina=0;
    IndexStr.fp_TP=full(sum(TPa_ind));
    IndexStr.fn_TP=full(sum(TPg_ind));
end

Dag_TP_temp=Dag_TP;
Dag_TP_temp(Dag_TP_temp==0)=inf;
Dag_TP_temp=min(Dag_TP_temp,[],2);
Dag_TP_temp(Dag_TP_temp==inf)=[]; 
if ~isempty(Dag_TP_temp)
    IndexStr.Dag_TP_mean=full(mean(Dag_TP_temp));
    IndexStr.Dag_TP_SD=full(std(Dag_TP_temp));
else
    IndexStr.Dag_TP_mean=NaN;
    IndexStr.Dag_TP_SD=NaN;
end

function IndexStr = CompareTraces(AM_G, r_G, R_G, AM_A, r_A, R_A, h_length,h_bp,h_tp)
% CompareTraces outputs parameters that compare the automated trace to the
%gold standard trace.
%set h_length=5l;h_bp=10;h_tp=10;
%R_G/A=zeros(size(r_G/A,1),1);

h_short_br=15;
h=max([h_length;h_bp;h_tp]);
s=100; % length threshold marking trees as corresponding

trace=CompareStructure(AM_G,r_G,R_G,AM_A,r_A,R_A,h,h_length,h_bp,h_tp,h_short_br);

AM_G(AM_G>0)=1;
AM_A(AM_A>0)=1;
AMlbl_trees_G=AM_G;
AMlbl_trees_A=AM_A;
%Delete lines above and uncomment lines below to allow for trees

% AMlbl_trees_G = LabelTreesAM(AM_G);
% AMlbl_trees_A = LabelTreesAM(AM_A);
Labels_G=unique(AMlbl_trees_G(AMlbl_trees_G>0));
Labels_A=unique(AMlbl_trees_A(AMlbl_trees_A>0));

% Parameters for individual traces
IndexStr.N_trees_G=length(Labels_G);
IndexStr.N_trees_A=length(Labels_A);
IndexStr.Lg=trace.Lg;
IndexStr.La=trace.La;
IndexStr.BPg=trace.BPg;
IndexStr.BPa=trace.BPa;
IndexStr.TPg=trace.TPg;
IndexStr.TPa=trace.TPa;

% Comparison of entire traces
IndexStr.corresponding_length_ing=trace.corresponding_length_ing;
IndexStr.corresponding_length_ina=trace.corresponding_length_ina;
IndexStr.fp_length=trace.fp_length;
IndexStr.fn_length=trace.fn_length;
IndexStr.Tortuosity_ag=IndexStr.corresponding_length_ina/IndexStr.corresponding_length_ing-1;
IndexStr.Dag_trace_full=trace.Dag_trace_full;
IndexStr.Dag_trace_mean=trace.Dag_trace_mean;
IndexStr.Dag_trace_SD=trace.Dag_trace_SD;

IndexStr.corresponding_BP_ing=trace.corresponding_BP_ing;
IndexStr.corresponding_BP_ina=trace.corresponding_BP_ina;
IndexStr.fp_BP=trace.fp_BP;
IndexStr.fn_BP=trace.fn_BP;
IndexStr.Dag_BP_full=trace.Dag_BP_full;
IndexStr.Dag_BP_mean=trace.Dag_BP_mean;
IndexStr.Dag_BP_SD=trace.Dag_BP_SD;

IndexStr.corresponding_TP_ing=trace.corresponding_TP_ing;
IndexStr.corresponding_TP_ina=trace.corresponding_TP_ina;
IndexStr.fp_TP=trace.fp_TP;
IndexStr.fn_TP=trace.fn_TP;
IndexStr.Dag_TP_full=trace.Dag_TP_full;
IndexStr.Dag_TP_mean=trace.Dag_TP_mean;
IndexStr.Dag_TP_SD=trace.Dag_TP_SD;

% Comparison of corresponding trees
IndexStr.Tree_Length_Recall=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);
IndexStr.Tree_BP_Recall=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);
IndexStr.Tree_TP_Recall=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);

IndexStr.Tree_Length_Precision=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);
IndexStr.Tree_BP_Precision=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);
IndexStr.Tree_TP_Precision=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);

IndexStr.Tree_Length_F_measure=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);
IndexStr.Tree_BP_F_measure=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);
IndexStr.Tree_TP_F_measure=nan(length(s),IndexStr.N_trees_G,IndexStr.N_trees_A);

NaN_temp=nan(IndexStr.N_trees_G,IndexStr.N_trees_A);
for k=1:length(s)
    TreesIndexStr=struct('Lg',NaN_temp,'La',NaN_temp,...
        'BPg',NaN_temp,'BPa',NaN_temp, ...
        'TPg',NaN_temp,'TPa',NaN_temp,...
        'corresponding_length_ing',NaN_temp,'corresponding_length_ina',NaN_temp, ...
        'fp_length',NaN_temp,'fn_length',NaN_temp,...
        'fn_BP',NaN_temp,'fp_BP',NaN_temp,...
        'fn_TP',NaN_temp,'fp_TP',NaN_temp, ...
        'Dag_trace_mean',NaN_temp,'Dag_trace_SD',NaN_temp,...
        'Dag_TP_mean',NaN_temp,'Dag_TP_SD',NaN_temp,...
        'Dag_BP_mean',NaN_temp,'Dag_BP_SD',NaN_temp, ...
        'Length_Recall',NaN_temp, 'Length_Precision',NaN_temp,...
        'BP_Recall',NaN_temp,'BP_Precision',NaN_temp,...
        'TP_Recall',NaN_temp,'TP_Precision',NaN_temp,...
        'Length_F_measure',NaN_temp,'BP_F_measure',NaN_temp,'TP_F_measure',NaN_temp);
    
    for i=1:length(Labels_G)
        [e1,e2]=find(AMlbl_trees_G==Labels_G(i));
        ee=unique([e1;e2]);
        AMg=AMlbl_trees_G(ee,ee);
        rg=r_G(ee,:);
        Rg=R_G(ee);
        for j=1:length(Labels_A),
            disp([i,j]);
            [e1,e2]=find(AMlbl_trees_A==Labels_A(j));
            ee=unique([e1;e2]);
            AMa=AMlbl_trees_A(ee,ee);
            ra=r_A(ee,:);
            Ra=R_A(ee);
            
            temp=CompareStructure(AMg, rg, Rg, AMa, ra, Ra, h,h_length,h_bp,h_tp,h_short_br);
            TreesIndexStr.Lg(i,j)=temp.Lg;
            TreesIndexStr.La(i,j)=temp.La;
            TreesIndexStr.BPg(i,j)=temp.BPg;
            TreesIndexStr.BPa(i,j)=temp.BPa;
            TreesIndexStr.TPg(i,j)=temp.TPg;
            TreesIndexStr.TPa(i,j)=temp.TPa;
            TreesIndexStr.Dag_trace_mean(i,j)=temp.Dag_trace_mean;
            TreesIndexStr.Dag_trace_SD(i,j)=temp.Dag_trace_SD;
            TreesIndexStr.corresponding_length_ing(i,j)=temp.corresponding_length_ing;
            TreesIndexStr.corresponding_length_ina(i,j)=temp.corresponding_length_ina;
            TreesIndexStr.corresponding_BP_ing(i,j)=temp.corresponding_BP_ing;
            TreesIndexStr.corresponding_BP_ina(i,j)=temp.corresponding_BP_ina;
            TreesIndexStr.corresponding_TP_ing(i,j)=temp.corresponding_TP_ing;
            TreesIndexStr.corresponding_TP_ina(i,j)=temp.corresponding_TP_ina;
            
            TreesIndexStr.fp_length(i,j)=temp.fp_length;
            TreesIndexStr.fn_length(i,j)=temp.fn_length;
            TreesIndexStr.fn_TP(i,j)=temp.fn_TP;
            TreesIndexStr.fp_TP(i,j)=temp.fp_TP;
            TreesIndexStr.Dag_TP_mean(i,j)=temp.Dag_TP_mean;
            TreesIndexStr.Dag_TP_SD(i,j)=temp.Dag_TP_SD;
            TreesIndexStr.fn_BP(i,j)=temp.fn_BP;
            TreesIndexStr.fp_BP(i,j)=temp.fp_BP;
            TreesIndexStr.Dag_BP_mean(i,j)=temp.Dag_BP_mean;
            TreesIndexStr.Dag_BP_SD(i,j)=temp.Dag_BP_SD;
            
            %Calculating the precision recall curve values for corresponding trees
            if (TreesIndexStr.corresponding_length_ing(i,j)+TreesIndexStr.corresponding_length_ina(i,j))/2>s(k)
                
                IndexStr.Tree_Length_Recall(k,i,j)=TreesIndexStr.corresponding_length_ing(i,j)./TreesIndexStr.Lg(i,j);
                IndexStr.Tree_BP_Recall(k,i,j)=TreesIndexStr.corresponding_BP_ing(i,j)./TreesIndexStr.BPg(i,j);
                IndexStr.Tree_TP_Recall(k,i,j)=TreesIndexStr.corresponding_TP_ing(i,j)./TreesIndexStr.TPg(i,j);
                
                IndexStr.Tree_Length_Precision(k,i,j)=TreesIndexStr.corresponding_length_ina(i,j)./TreesIndexStr.La(i,j);
                IndexStr.Tree_BP_Precision(k,i,j)=TreesIndexStr.corresponding_BP_ina(i,j)./TreesIndexStr.BPa(i,j);
                IndexStr.Tree_TP_Precision(k,i,j)=TreesIndexStr.corresponding_TP_ina(i,j)./TreesIndexStr.TPa(i,j);
                
                IndexStr.Tree_Length_F_measure(k,i,j)=2./(1./TreesIndexStr.Length_Recall(i,j)+1./TreesIndexStr.Length_Precision(i,j));
                IndexStr.Tree_BP_F_measure(k,i,j)=2./(1./TreesIndexStr.BP_Recall(i,j)+1./TreesIndexStr.BP_Precision(i,j));
                IndexStr.Tree_TP_F_measure(k,i,j)=2./(1./TreesIndexStr.TP_Recall(i,j)+1./TreesIndexStr.TP_Precision(i,j));
            end
        end
    end
end
end

%CompareStructure compares 2 given traces, trees or branches depending on
%input matrices.
function IndexStr = CompareStructure(AMg, rg, Rg, AMa, ra, Ra, h, h_length ,h_bp, h_tp, h_short_br)

%removing AM labels
AMg(AMg>0)=1;
AMa(AMa>0)=1;

%Calculating total length of the structures
[ig,jg]=find(AMg);
lllg=AMg;
lllg(AMg>0)=sum((rg(ig,:)-rg(jg,:)).^2,2).^0.5;
lg=sum(lllg,1)./2;
IndexStr.Lg=full(sum(lg));

[ia,ja]=find(AMa);
llla=AMa; llla(AMa>0)=sum((ra(ia,:)-ra(ja,:)).^2,2).^0.5;
la=sum(llla,2)./2;
IndexStr.La=full(sum(la));

%Finding all Branch Points in the structure
BPg_ind=(sum(AMg)>2);
IndexStr.BPg=full(sum(BPg_ind));
BPa_ind=(sum(AMa)>2);
IndexStr.BPa=full(sum(BPa_ind));

%Finding Terminal Points in the structure
TPg_ind=(sum(AMg)==1);
IndexStr.TPg=full(sum(TPg_ind));
TPa_ind=(sum(AMa)==1);
IndexStr.TPa=full(sum(TPa_ind));

%Calculation of the node to node separation matrix for the structure

temp=(size(AMg,1)+size(AMa,1))*100*ceil(2*h);
Dag_ijv=zeros(temp,3);
Dag_minsparseinds=zeros(size(AMg,1),2);
count=1;
for j=1:size(AMg,1)
    d=((ra(:,1)-rg(j,1)).^2+(ra(:,2)-rg(j,2)).^2+(ra(:,3)-rg(j,3)).^2).^0.5; %Euclidean distance between vertices
    ind=find(d<=h);
    d=d(ind);
    d(d==0)=0.001; % Fixed to avoid having to deal Storage of zeros in sparse format
    if ~isempty(ind)
        countend=count+length(ind)-1;
        Dag_ijv(count:countend,:)=[j.*ones(length(ind),1),ind,d];
        [d_min,ind_min]=min(Dag_ijv(count:countend,3));
        ind_min=Dag_ijv(ind_min+count-1,2);
        Dag_minsparseinds(j,:)=[ind_min,d_min];
        count=countend+1;
    end
end
Dag_ijv=Dag_ijv(1:count-1,:);
temp=Dag_minsparseinds(:,2)>h_length; %finding points that lie within h_length
Dag_minsparseinds(temp,2)=0; %sparse indexing of Dag_min will ignore these
temp=(Dag_minsparseinds(:,1)==0); %finding points that were not assigned values
Dag_minsparseinds(temp,1)=1; %Arbitrary value which will be ignored while sparse indexing

j=1:size(AMg,1);
Dag_min=sparse(j,Dag_minsparseinds(j,1),Dag_minsparseinds(j,2),size(AMg,1),size(AMa,1));
Dag=sparse(Dag_ijv(:,1),Dag_ijv(:,2),Dag_ijv(:,3),size(AMg,1),size(AMa,1));

% Find lengths associated with the Terminal and Branch Points
AMlblg=LabelBranchesAM(AMg);
BranchLengthsg=BranchLengthsAM(AMlblg,rg);
TPg_length=BranchLengthsg(sum(AMlblg(:,TPg_ind)));

temp=AMlblg(:,BPg_ind);
temp2=sort(nonzeros(temp));
interbranches=unique(temp2(~([0;temp2(1:end-1)]-temp2)));
BranchLengthsg(interbranches)=inf;
temp_lengths=inf(size(temp));
temp_lengths(temp>0)=BranchLengthsg(temp(temp>0));
BPg_length=min(temp_lengths);

AMlbla=LabelBranchesAM(AMa);
BranchLengthsa=BranchLengthsAM(AMlbla,ra);
TPa_length=BranchLengthsa(sum(AMlbla(:,TPa_ind)));

temp=AMlbla(:,BPa_ind);
temp2=sort(nonzeros(temp));
interbranches=unique(temp2(~([0;temp2(1:end-1)]-temp2)));
BranchLengthsa(interbranches)=inf;
temp_lengths=inf(size(temp));
temp_lengths(temp>0)=BranchLengthsa(temp(temp>0));
BPa_length=min(temp_lengths);

NaNg=zeros(size(Dag,1),1);
NaNg(TPg_ind)=(TPg_length<h_short_br);
NaNg(BPg_ind)=(BPg_length<h_short_br);
NaNg(NaNg==1)=NaN;

NaNa=zeros(1,size(Dag,2));
NaNa(TPa_ind)=(TPa_length<h_short_br);
NaNa(BPa_ind)=(BPa_length<h_short_br);
NaNa(NaNa==1)=NaN;

NaNag=sparse(NaNg*NaNa);

%Calculation of distance between the structures, and lengths of corresponding,
%false positive/negative parts of the structure
IndexStr.Dag_trace_full=nonzeros(Dag_min(:));
IndexStr.Dag_trace_mean=full((lg*Dag_min*la)/(lg*(Dag_min>0)*la));
IndexStr.Dag_trace_SD=full(((lg*Dag_min.^2*la)/(lg*(Dag_min>0)*la)-IndexStr.Dag_trace_mean^2)^0.5);
IndexStr.corresponding_length_ing=full(lg*(sum(Dag,2)>0));
IndexStr.corresponding_length_ina=full((sum(Dag,1)>0)*la);
IndexStr.fp_length=full((sum(Dag,1)==0)*la);
IndexStr.fn_length=full(lg*(sum(Dag,2)==0));

% Calculation of distance between corresponding branch points and number of
% corresponding, false positive/negative branch points
Dag_BP=Dag(BPg_ind,BPa_ind);
Dag_BP(Dag_BP>h_bp)=0;
NaNag_BP=NaNag(BPg_ind,BPa_ind);
Dag_BP(Dag_BP==0)=NaNag_BP(Dag_BP==0);

if ~isempty(Dag_BP)
    IndexStr.corresponding_BP_ing=full(sum(nansum(Dag_BP,2)>0));
    IndexStr.corresponding_BP_ina=full(sum(nansum(Dag_BP,1)>0));
    IndexStr.fp_BP=full(sum(nansum(Dag_BP,1)==0))-sum(prod((double(isnan(Dag_BP))),1));
    IndexStr.fn_BP=full(sum(nansum(Dag_BP,2)==0))-sum(prod((double(isnan(Dag_BP))),2));
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
    IndexStr.Dag_BP_full=Dag_BP_temp(~isnan(Dag_BP_temp));
    IndexStr.Dag_BP_mean=full(nanmean(Dag_BP_temp));
    IndexStr.Dag_BP_SD=full(nanstd(Dag_BP_temp));
else
    IndexStr.Dag_BP_full=NaN;
    IndexStr.Dag_BP_mean=NaN;
    IndexStr.Dag_BP_SD=NaN;
end

% Calculation of distance between corresponding terminal points and number of
% corresponding, false positive/negative terminal points
Dag_TP=Dag(TPg_ind,TPa_ind);
Dag_TP(Dag_TP>h_tp)=0;
NaNag_TP=NaNag(TPg_ind,TPa_ind);
Dag_TP(Dag_TP==0)=NaNag_TP(Dag_TP==0);

if ~isempty(Dag_TP)
    IndexStr.corresponding_TP_ing=full(sum(nansum(Dag_TP,2)>0));
    IndexStr.corresponding_TP_ina=full(sum(nansum(Dag_TP,1)>0));
    IndexStr.fp_TP=full(sum(nansum(Dag_TP,1)==0))-sum(prod((double(isnan(Dag_TP))),1));
    IndexStr.fn_TP=full(sum(nansum(Dag_TP,2)==0))-sum(prod((double(isnan(Dag_TP))),2));
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
    IndexStr.Dag_TP_full=Dag_TP_temp(~isnan(Dag_TP_temp));
    IndexStr.Dag_TP_mean=full(nanmean(Dag_TP_temp));
    IndexStr.Dag_TP_SD=full(nanstd(Dag_TP_temp));
else
    IndexStr.Dag_TP_mean=NaN;
    IndexStr.Dag_TP_SD=NaN;
end

end

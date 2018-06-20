function Out = CompareTraces_byLength(AMg, rg, AMa, ra, h)
%Optimized to disregard trees, and calculate only inter-trace distance and
%length based measures

%removing AM labels
AMg(AMg>0)=1;
AMa(AMa>0)=1;

%Calculating total length of the structures
[ig,jg]=find(AMg);
lllg=AMg;
lllg(AMg>0)=sum((rg(ig,:)-rg(jg,:)).^2,2).^0.5;
lg=sum(lllg,1)./2;
lg=full(lg);
Out.Lg=sum(lg);

[ia,ja]=find(AMa);
llla=AMa; llla(AMa>0)=sum((ra(ia,:)-ra(ja,:)).^2,2).^0.5;
la=sum(llla,2)./2;
la=full(la);
Out.La=sum(la);

%All to all inter-vertex distances:
Dag=(bsxfun(@minus,ra(:,1),rg(:,1)').^2+bsxfun(@minus,ra(:,2),rg(:,2)').^2+bsxfun(@minus,ra(:,3),rg(:,3)').^2).^0.5;
 Dag(Dag>h)=nan;%Distances exceeding threshold set to nan

[Out.Dag_trace_full,Out.bestmatch_to_g_in_a]=nanmin(Dag,[],1);
Out.bestmatch_to_g_in_a(isnan(Out.Dag_trace_full))=nan;
[temp,Out.bestmatch_to_a_in_g]=nanmin(Dag,[],2);
Out.bestmatch_to_a_in_g(isnan(temp))=nan;
Out.Dag_trace_full(isnan(Out.Dag_trace_full))=[];

%Indices of trace vertices that are further than threshold from any
%vertices in the other trace:
far_in_g=isnan(nanmin(Dag,[],1));
far_in_a=isnan(nanmin(Dag,[],2));

%Calculation of distance between the structures, and lengths of corresponding,
%false positive/negative parts
Out.corresponding_length_ing=lg*~far_in_g(:);
Out.corresponding_length_ina=~far_in_a(:)'*la;
Out.fp_length=far_in_a(:)'*la;
Out.fn_length=lg*far_in_g(:);
end

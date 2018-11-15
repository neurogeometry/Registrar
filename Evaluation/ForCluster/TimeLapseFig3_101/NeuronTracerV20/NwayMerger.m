% This function determines the costs of different n-way merger scenarios.
% Matrices Dist,Cos,Offset,MeanI,CVI contain cost components for pair-wise mergers 
% Cell array MergerAM contains sparse AMs for different merging scenarios.
% The corresponding Costs are ordered from lowest to highest.

function [Costs MergerAM]=NwayMerger(TipLabels,Dist,Overrun,Offset,Cos,Parameters)

d=Parameters.Classification.DistanceWeight;
ov=Parameters.Classification.OverrunWeight;
of=Parameters.Classification.OffsetWeight;
c=Parameters.Classification.AnglesWeight;
g=Parameters.Classification.TerminalsWeight;

Cos2=-1; Cos3=-0.5; Cos4=0;

N=size(Dist,1); % number of points in the cluster to be merged

Bins={};
if N==1
    Bins{1}=1;
elseif N==2
    Bins{1}=[1,1];
    Bins{2}=2;
elseif N==3
    Bins{1}=[1,1,1];
    Bins{2}=[2,1];
    Bins{3}=3;
elseif N==4
    Bins{1}=[1,1,1,1];
    Bins{2}=[2,1,1];
    Bins{3}=[2,2];
    Bins{4}=[3,1];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{5}=[4];
    end
elseif N==5
    Bins{1}=[1,1,1,1,1];
    Bins{2}=[2,1,1,1];
    Bins{3}=[2,2,1];
    Bins{4}=[3,1,1];
    Bins{5}=[3,2];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{6}=[4,1];
    end
elseif N==6
    Bins{1}=[1,1,1,1,1,1];
    Bins{2}=[2,1,1,1,1];
    Bins{3}=[2,2,1,1];
    Bins{4}=[2,2,2];
    Bins{5}=[3,1,1,1];
    Bins{6}=[3,2,1];
    Bins{7}=[3,3];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{8}=[4,1,1];
        Bins{9}=[4,2];
    end
elseif N==7
    Bins{1}=[1,1,1,1,1,1,1];
    Bins{2}=[2,1,1,1,1,1];
    Bins{3}=[2,2,1,1,1];
    Bins{4}=[2,2,2,1];
    Bins{5}=[3,1,1,1,1];
    Bins{6}=[3,2,1,1];
    Bins{7}=[3,2,2];
    Bins{8}=[3,3,1];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{9}=[4,1,1,1];
        Bins{10}=[4,2,1];
        Bins{11}=[4,3];
    end
elseif N==8
    Bins{1}=[1,1,1,1,1,1,1,1];
    Bins{2}=[2,1,1,1,1,1,1];
    Bins{3}=[2,2,1,1,1,1];
    Bins{4}=[2,2,2,1,1];
    Bins{5}=[2,2,2,2];
    Bins{6}=[3,1,1,1,1,1];
    Bins{7}=[3,2,1,1,1];
    Bins{8}=[3,2,2,1];
    Bins{9}=[3,3,1,1];
    Bins{10}=[3,3,2];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{11}=[4,1,1,1,1];
        Bins{12}=[4,2,1,1];
        Bins{13}=[4,2,2];
        Bins{14}=[4,3,1];
        Bins{15}=[4,4];
    end
elseif N==9
    Bins{1}=[1,1,1,1,1,1,1,1,1];
    Bins{2}=[2,1,1,1,1,1,1,1];
    Bins{3}=[2,2,1,1,1,1,1];
    Bins{4}=[2,2,2,1,1,1];
    Bins{5}=[2,2,2,2,1];
    Bins{6}=[3,1,1,1,1,1,1];
    Bins{7}=[3,2,1,1,1,1];
    Bins{8}=[3,2,2,1,1];
    Bins{9}=[3,2,2,2];
    Bins{10}=[3,3,1,1,1];
    Bins{11}=[3,3,2,1];
    Bins{12}=[3,3,3];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{13}=[4,1,1,1,1,1];
        Bins{14}=[4,2,1,1,1];
        Bins{15}=[4,2,2,1];
        Bins{16}=[4,3,1,1];
        Bins{17}=[4,3,2];
        Bins{18}=[4,4,1];
    end
elseif N==10
    Bins{1}=[1,1,1,1,1,1,1,1,1,1];
    Bins{2}=[2,1,1,1,1,1,1,1,1];
    Bins{3}=[2,2,1,1,1,1,1,1];
    Bins{4}=[2,2,2,1,1,1,1];
    Bins{5}=[2,2,2,2,1,1];
    Bins{6}=[2,2,2,2,2];
    Bins{7}=[3,1,1,1,1,1,1,1];
    Bins{8}=[3,2,1,1,1,1,1];
    Bins{9}=[3,2,2,1,1,1];
    Bins{10}=[3,2,2,2,1];
    Bins{11}=[3,3,1,1,1,1];
    Bins{12}=[3,3,2,1,1];
    Bins{13}=[3,3,2,2];
    Bins{14}=[3,3,3,1];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{15}=[4,1,1,1,1,1,1];
        Bins{16}=[4,2,1,1,1,1];
        Bins{17}=[4,2,2,1,1];
        Bins{18}=[4,2,2,2];
        Bins{19}=[4,3,1,1,1];
        Bins{20}=[4,3,2,1];
        Bins{21}=[4,3,3];
        Bins{22}=[4,4,1,1];
        Bins{23}=[4,4,2];
    end
elseif N==11
    Bins{1}=[1,1,1,1,1,1,1,1,1,1,1];
    Bins{2}=[2,1,1,1,1,1,1,1,1,1];
    Bins{3}=[2,2,1,1,1,1,1,1,1];
    Bins{4}=[2,2,2,1,1,1,1,1];
    Bins{5}=[2,2,2,2,1,1,1];
    Bins{6}=[2,2,2,2,2,1];   
    Bins{7}=[3,1,1,1,1,1,1,1,1];
    Bins{8}=[3,2,1,1,1,1,1,1];
    Bins{9}=[3,2,2,1,1,1,1];
    Bins{10}=[3,2,2,2,1,1];
    Bins{11}=[3,2,2,2,2];
    Bins{12}=[3,3,1,1,1,1,1];
    Bins{13}=[3,3,2,1,1,1];
    Bins{14}=[3,3,2,2,1];
    Bins{15}=[3,3,3,1,1];
    Bins{16}=[3,3,3,2];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{17}=[4,1,1,1,1,1,1,1];
        Bins{18}=[4,2,1,1,1,1,1];
        Bins{19}=[4,2,2,1,1,1];
        Bins{20}=[4,2,2,2,1];
        Bins{21}=[4,3,1,1,1,1];
        Bins{22}=[4,3,2,1,1];
        Bins{23}=[4,3,2,2];
        Bins{24}=[4,3,3,1];
        Bins{25}=[4,4,1,1,1];
        Bins{26}=[4,4,2,1];
        Bins{27}=[4,4,3];
    end
elseif N==12
    Bins{1}=[1,1,1,1,1,1,1,1,1,1,1,1];
    Bins{2}=[2,1,1,1,1,1,1,1,1,1,1];
    Bins{3}=[2,2,1,1,1,1,1,1,1,1];
    Bins{4}=[2,2,2,1,1,1,1,1,1];
    Bins{5}=[2,2,2,2,1,1,1,1];
    Bins{6}=[2,2,2,2,2,1,1];
    Bins{7}=[2,2,2,2,2,2];
    Bins{8}=[3,1,1,1,1,1,1,1,1,1];
    Bins{9}=[3,2,1,1,1,1,1,1,1];
    Bins{10}=[3,2,2,1,1,1,1,1];
    Bins{11}=[3,2,2,2,1,1,1];
    Bins{12}=[3,2,2,2,2,1];
    Bins{13}=[3,3,1,1,1,1,1,1];
    Bins{14}=[3,3,2,1,1,1,1];
    Bins{15}=[3,3,2,2,1,1];
    Bins{16}=[3,3,2,2,2];
    Bins{17}=[3,3,3,1,1,1];
    Bins{18}=[3,3,3,2,1];
    Bins{19}=[3,3,3,3];
    if Parameters.BranchMerger.FourWayMerger==1
        Bins{20}=[4,1,1,1,1,1,1,1,1];
        Bins{21}=[4,2,1,1,1,1,1,1];
        Bins{22}=[4,2,2,1,1,1,1];
        Bins{23}=[4,2,2,2,1,1];
        Bins{24}=[4,2,2,2,2];
        Bins{25}=[4,3,1,1,1,1,1];
        Bins{26}=[4,3,2,1,1,1];
        Bins{27}=[4,3,2,2,1];
        Bins{28}=[4,3,3,1,1];
        Bins{29}=[4,3,3,2];
        Bins{30}=[4,4,1,1,1,1];
        Bins{31}=[4,4,2,1,1];
        Bins{32}=[4,4,2,2];
        Bins{33}=[4,4,3,1];
        Bins{34}=[4,4,4];
    end
end

MergerAM=[];
Costs=[];

% Merging constraints
Forbid=(Dist>Parameters.BranchMerger.MaxDistance | (TipLabels'*ones(1,length(TipLabels))-ones(length(TipLabels),1)*TipLabels)==0); %| Overrun>5 
Force=[];
% [temp1,temp2]=find(Dist==0);
% Force=unique([temp1,temp2]);

for i=1:length(Bins)
    A=permutations(Bins{i},Forbid,Force);
    if ~isempty(A)
        MergerAMtemp=zeros(N,N,size(A,1));
        Costs_temp=zeros(size(A,1),1);
        for j=1:length(Bins{i})
            if Bins{i}(j)==2
                temp_ij=sort(A(:,sum( Bins{i}(1:j))-1:sum( Bins{i}(1:j))),2);
                MergerAMtemp(sub2ind(size(MergerAMtemp),temp_ij(:,1),temp_ij(:,2),(1:size(A,1))'))=1;
                
                temp_ind=sub2ind([N,N],temp_ij(:,1),temp_ij(:,2));
                Costs_temp=Costs_temp+d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos2);
                                
            elseif Bins{i}(j)==3
                temp_ij=sort(A(:,sum( Bins{i}(1:j))-2:sum( Bins{i}(1:j))),2);
                MergerAMtemp(sub2ind(size(MergerAMtemp),temp_ij(:,1),temp_ij(:,2),(1:size(A,1))'))=1;
                MergerAMtemp(sub2ind(size(MergerAMtemp),temp_ij(:,1),temp_ij(:,3),(1:size(A,1))'))=1;
                % !!!!! triangle is not closed to avoid forming a loop
                
                temp_ind=sub2ind([N,N],temp_ij(:,1),temp_ij(:,2));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos3))./3;
                temp_ind=sub2ind([N,N],temp_ij(:,1),temp_ij(:,3));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos3))./3;
                temp_ind=sub2ind([N,N],temp_ij(:,2),temp_ij(:,3));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos3))./3;
                                
            elseif Bins{i}(j)==4
                temp_ij=sort(A(:,sum( Bins{i}(1:j))-3:sum( Bins{i}(1:j))),2);
                MergerAMtemp(sub2ind(size(MergerAMtemp),temp_ij(:,1),temp_ij(:,2),(1:size(A,1))'))=1;
                MergerAMtemp(sub2ind(size(MergerAMtemp),temp_ij(:,1),temp_ij(:,3),(1:size(A,1))'))=1;
                MergerAMtemp(sub2ind(size(MergerAMtemp),temp_ij(:,1),temp_ij(:,4),(1:size(A,1))'))=1;
                % !!!!! square is not closed to avoid forming a loop
                
                temp_ind=sub2ind([N,N],temp_ij(:,1),temp_ij(:,2));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos4))./6;
                temp_ind=sub2ind([N,N],temp_ij(:,1),temp_ij(:,3));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos4))./6;
                temp_ind=sub2ind([N,N],temp_ij(:,1),temp_ij(:,4));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos4))./6;
                temp_ind=sub2ind([N,N],temp_ij(:,2),temp_ij(:,3));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos4))./6;
                temp_ind=sub2ind([N,N],temp_ij(:,2),temp_ij(:,4));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos4))./6;
                temp_ind=sub2ind([N,N],temp_ij(:,3),temp_ij(:,4));
                Costs_temp=Costs_temp+(d.*Dist(temp_ind)+ov.*Overrun(temp_ind)+of.*Offset(temp_ind)+c.*abs(Cos(temp_ind)-Cos4))./6;
            end
        end
        MergerAM=cat(3,MergerAM,MergerAMtemp);
        Costs=[Costs;Costs_temp+g.*nnz(Bins{i}==1)];
    end
end

% Sort the mergers by cost
[Costs,Cost_ind]=sort(Costs);
MergerAM=MergerAM(:,:,Cost_ind);

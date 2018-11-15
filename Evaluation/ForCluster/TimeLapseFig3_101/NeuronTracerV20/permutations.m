% This function generates all different permutations of placing distinguishable 
% objects into indistinguishable bins of different sizes. Matrix Forbid contains 
% information about pairs of objects which should not be placed in the same bin. 
% Vector Force contains objects which must not be placed in bins of size 1.

function A=permutations(BinSizes,Forbid,Force)

Nbins=length(BinSizes);
Nballs=sum(BinSizes);

if isempty(Forbid)
    Forbid=zeros(Nballs,Nballs);
end

A=zeros(1,Nballs);
R=rand(1,Nballs);
i=0;
while i<Nballs && ~isempty(A)
    i=i+1;
    A=add_a_ball(A,BinSizes,i,Forbid,Force);
    
    Rtemp=A;
    Rtemp(A==0)=rand;
    Rtemp(A~=0)=R(A(A~=0));
    
    B=zeros(size(A,1),1);
    start=1;
    for j=1:Nbins
        B=B+prod(Rtemp(:,start:start+BinSizes(j)-1),2);
        start=start+BinSizes(j);
    end
    [~, m, ~]=unique(fix(B*10^10));
    A=A(m,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A_new=add_a_ball(A,BinSizes,b,Forbid,Force)

Nbins=length(BinSizes);
Nballs=sum(BinSizes);
begs=cumsum(BinSizes)-BinSizes+1;
ends=cumsum(BinSizes);

if ismember(b,Force)
    temp=zeros(size(A));
    temp(:,begs(begs==ends))=1;
    A(temp & A==0)=-1;
end

B=zeros(size(A,1),1);
for i=1:Nbins
    B=B+(prod(A(:,begs(i):ends(i)),2)==0);
end

A_new=zeros(sum(B),Nballs);
count0=1;
for j=1:size(A,1)
    count=0;
    A_temp=zeros(B(j),Nballs);
    for i=1:Nbins
        temp=A(j,begs(i):ends(i));
        ind=find(temp==0,1,'first');
        if ~isempty(ind) 
            if sum(Forbid(b,nonzeros(temp)))==0
                count=count+1;
                A_temp(count,:)=A(j,:);
                A_temp(count,begs(i)+ind-1)=b;
            end
        end
    end
    A_temp=A_temp(1:count,:);
    A_new(count0:count0-1+size(A_temp,1),:)=A_temp;
    count0=count0+size(A_temp,1);
end
A_new=A_new(1:count0-1,:);
A_new(A_new==-1)=0;

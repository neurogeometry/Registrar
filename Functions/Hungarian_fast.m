% Implementation of Hungarian algorithm according to Munkres 1957 and  Burgeois 1971
% Cost is the cost matrix. It may be rectangular and may not contain Infs or NaNs
% UnassignedWorkerCost and UnassignedJobCost are vectors (scalars or empty) specifying the
% costs of not assigning workers and jobs respectively. They may not contain Infs or NaNs.
% If UnassignedWorkerCost == [], UnassignedJobCost must == [] as well, in which case the 
% maximum possible number of assignments must be made.
% AM is the assignment matrix and E is the assignment cost

% Steps of the algorithm are based on http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
% and F. Burgeois and J.-C. Lasalle. An extension of the Munkres algorithm for the assignment problem to rectangular matrices. Communications of the ACM, 142302-806, 1971.

% Step 1:   For each row of the matrix, find the smallest element and subtract it from every element in its row.  Go to Step 2.
% Step 2:   Find a zero (Z) in the resulting matrix.  If there is no starred zero in its row or column, star Z. 
%           Repeat for each element in the matrix. Go to Step 3.
% Step 3:   Cover each column containing a starred zero. If K columns are covered, the starred zeros describe a complete set of unique assignments.
%           In this case, Go to DONE, otherwise, Go to Step 4.
% Step 4:   Find a noncovered zero and prime it. If there is no starred zero in the row containing this primed zero, Go to Step 5.
%           Otherwise, cover this row and uncover the column containing the starred zero. Continue in this manner (i.e. Find a 
%           noncovered zero and prime it...) until there are no uncovered zeros left. Save the smallest uncovered value and Go to Step 6.
% Step 5:   Construct a series of alternating primed and starred zeros as follows.  Let Z0 represent the uncovered primed zero found in Step 4.
%           Let Z1 denote the starred zero in the column of Z0 (if any). Let Z2 denote the primed zero in the row of Z1 (there will always be one).
%           Continue until the series terminates at a primed zero that has no starred zero in its column.  Unstar each starred zero of the series, 
%           star each primed zero of the series, erase all primes and uncover every line in the matrix.  Return to Step 3.
% Step 6:   Add the value found in Step 4 to every element of each covered row, and subtract it from every element of each uncovered column.  
%           Return to Step 4 without altering any stars, primes, or covered lines.

function [AM,E]=Hungarian_fast(Cost,UnassignedWorkerCost,UnassignedJobCost)

validateattributes(Cost,{'numeric'},{'real','nonnan','nonempty','finite','2d'},'assignDetectionsToTracks','Cost');
if nargin==3
    validateattributes(UnassignedWorkerCost,{'numeric'},{'real','nonnan','finite'},'assignDetectionsToTracks','UnassignedWorkerCost');
    validateattributes(UnassignedJobCost,{'numeric'},{'real','nonnan','finite'},'assignDetectionsToTracks','UnassignedJobCost');
elseif nargin==1
    UnassignedWorkerCost=[];
    UnassignedJobCost=[];
end
assert(isempty(UnassignedWorkerCost) || isscalar(UnassignedWorkerCost) || all(size(UnassignedWorkerCost)==size(Cost(:,1))),'UnassignedWorkerCost must be empty, scalar, or vector of size(Cost(:,1))')
assert(isempty(UnassignedJobCost) || isscalar(UnassignedJobCost) || all(size(UnassignedJobCost)==size(Cost(1,:))),'UnassignedJobCost must be empty, scalar, or vector of size(Cost(1,:))')
assert((~isempty(UnassignedJobCost) && ~isempty(UnassignedWorkerCost)) || (isempty(UnassignedJobCost) && isempty(UnassignedWorkerCost)),'If UnassignedJobCost is empty, UnassignedWorkerCost must be empty as well')

if isscalar(UnassignedWorkerCost)
    UnassignedWorkerCost=UnassignedWorkerCost.*ones(size(Cost,1),1);
end
if isscalar(UnassignedJobCost)
    UnassignedJobCost=UnassignedJobCost.*ones(1,size(Cost,2));
end   

% works faster if size(Cost,1)>=size(Cost,2) 
Cost_transpose=false;
if size(Cost,1)<size(Cost,2)
    Cost=Cost';
    temp=UnassignedWorkerCost;
    UnassignedWorkerCost=UnassignedJobCost';
    UnassignedJobCost=temp';
    Cost_transpose=true;
end

% C is padded cost
if isempty(UnassignedWorkerCost) && isempty(UnassignedJobCost)
    C=Cost; 
elseif ~isempty(UnassignedWorkerCost) && ~isempty(UnassignedJobCost)
    BigNumber=max([max(Cost(:)),max(UnassignedWorkerCost),max(UnassignedJobCost)])*(sum(size(Cost)))^2+1;
    %C=[Cost,BigNumber.*ones(size(Cost,1),size(Cost,1))+diag(UnassignedWorkerCost-BigNumber);BigNumber.*ones(size(Cost,2),size(Cost,2))+diag(UnassignedJobCost-BigNumber),eps.*rand(size(Cost,2),size(Cost,1))];

    % Fast alternative
    SmallNumber=eps/(sum(size(Cost)))^2/10;
    C=[Cost,BigNumber.*ones(size(Cost,1),size(Cost,1))+diag(UnassignedWorkerCost-BigNumber);BigNumber.*ones(size(Cost,2),size(Cost,2))+diag(UnassignedJobCost-BigNumber),SmallNumber.*rand(size(Cost,2),size(Cost,1))];
end

N=size(C);
K=min(N);

% Step 1
if N(1)>N(2)
    C=C-ones(N(1),1)*min(C,[],1);
elseif N(1)<N(2)
    C=C-min(C,[],2)*ones(1,N(2));
else
    C=C-ones(N(1),1)*min(C,[],1);
    C=C-min(C,[],2)*ones(1,N(2));
end

% Step 2
Z=sparse(C==0); % matrix of all zeros
S=logical(sparse(N(1),N(2))); % matrix of starred zeros
P=logical(sparse(N(1),N(2))); % matrix of primed zeros

[i,j]=find(C==0);
while ~isempty(i)
    S(i(1),j(1))=true; % starred zero
    remind=(i==i(1) | j==j(1));
    i(remind)=[];
    j(remind)=[];
end

% Step 3
hlines=false(N(1),1);
vlines=false(1,N(2));
vlines(any(S,1))=true;
while sum(vlines)<K
    % Step 4
    [i_nonc,j_nonc]=find(Z(~hlines,~vlines),1); % noncovered zero
    if ~isempty(i_nonc)
        ind=find(hlines==0,i_nonc);
        i_nonc=ind(end);
        ind=find(vlines==0,j_nonc);
        j_nonc=ind(end);
    end
    
    if ~isempty(i_nonc)
        P(i_nonc,j_nonc)=true; % primed zero
        temp_ind=find(S(i_nonc,:),1);
        if isempty(temp_ind)
            % Step 5
            Z0=[i_nonc,j_nonc]; % primed zero
            inds=Z0(1)+N(1)*(Z0(2)-1);
            if ~isempty(find(S(:,Z0(2)),1))
                Z1=[find(S(:,Z0(2)),1),Z0(2)]; % starred zero
                inds=[inds,Z1(1)+N(1)*(Z0(2)-1)];
            else
                Z1=[];
            end
            while ~isempty(Z1)
                Z0=[Z1(1),find(P(Z1(1),:),1)]; % primed zero
                inds=[inds,Z0(1)+N(1)*(Z0(2)-1)];
                if ~isempty(find(S(:,Z0(2)),1))
                    Z1=[find(S(:,Z0(2)),1),Z0(2)]; % starred zero
                    inds=[inds,Z1(1)+N(1)*(Z0(2)-1)];
                else
                    Z1=[];
                end
            end
            if length(inds)>1
                S(inds(2:2:end-1))=false;
            end
            S(inds(1:2:end))=true;
            P=logical(logical(sparse(N(1),N(2))));
            hlines=false(N(1),1);
            vlines(any(S,1))=true;
        else
            hlines(i_nonc)=true;
            vlines(temp_ind)=false;
        end
    else
        % Step 6
        temp1=C(~hlines,~vlines);
        [M,ind]=min(temp1(:));
        [i_temp,j_temp]=ind2sub(size(temp1),ind);
        ind=find(hlines==0,i_temp);
        i_temp=ind(end);
        ind=find(vlines==0,j_temp);
        j_temp=ind(end);
        Z(i_temp,j_temp)=true;

        C(hlines,vlines)=C(hlines,vlines)+M;
        C(~hlines,~vlines)=temp1-M;
        
        Z(hlines,vlines)=false;
        S(hlines,vlines)=false;
        P(hlines,vlines)=false;
    end
end

AM=S(1:size(Cost,1),1:size(Cost,2));
E=sum(Cost(AM));
if ~isempty(UnassignedWorkerCost)
    E=E+sum(UnassignedWorkerCost(~any(AM,2)))+sum(UnassignedJobCost(~any(AM,1)));
end

if Cost_transpose
    AM=AM';
end


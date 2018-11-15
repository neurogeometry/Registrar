function [AM rt I_snake]=Snake_subd_AM_Sparse(Im,sizeIm,AMlbl,r)

pointsperum=1;%/2;
Nstep=40;%/4;
alpha=0.1;%/4;
betta=0.1;%*4; %1/8/alpha/8;
bpgrad=0;

rt=r';
AM = spones(AMlbl+AMlbl');

N=length(AM);
[it,jt]=find(AM);
lll=AM; lll(AM>0)=sum((rt(:,it)-rt(:,jt)).^2).^0.5;

format short g
full([N, sum(lll(:))/2])

Im=Im./max(Im(:));
sigmax=1;
sigmay=1;
sigmaz=1;

sigx=1.0; %0.5
sigy=1.0; %0.5
sigz=1.0; %0.5
deltax=fix(2*sigx);
deltay=fix(2*sigy);
deltaz=fix(2*sigz);
S=[2*deltax,2*deltay,2*deltaz]+1;
[xtemp,ytemp,ztemp]=ind2sub(S,1:prod(S));

for f1=1:Nstep
    % resegment
    [AM rt] = Divide_Segments_AS(AM,rt,pointsperum);
    [AM rt] = Merge_Segments_AS(AM,rt,pointsperum);
    [it,jt]=find(AM);
    N=length(AM);
    lll=AM; lll(AM>0)=sum((rt(:,it)-rt(:,jt)).^2).^0.5;
    notendp=(sum(AM)>1);
    bp=(sum(AM)>2);
    BP=zeros(1,N); BP(bp)=1; BP=[1;1;1]*BP;
    
    Fxtemp=zeros(prod(S),N); Fytemp=zeros(prod(S),N); Fztemp=zeros(prod(S),N);
    Itemp=zeros(prod(S),N); Dg=zeros(3,N);

    Gx=exp(-((rt(1,:)-fix(rt(1,:)))'*ones(1,prod(S))+deltax-ones(N,1)*xtemp).^2./2./sigx^2);
    Gy=exp(-((rt(2,:)-fix(rt(2,:)))'*ones(1,prod(S))+deltay-ones(N,1)*ytemp).^2./2./sigy^2);
    Gz=exp(-((rt(3,:)-fix(rt(3,:)))'*ones(1,prod(S))+deltaz-ones(N,1)*ztemp).^2./2./sigz^2);
    G=Gx.*Gy.*Gz;
    G=G./(sum(G,2)*ones(1,prod(S)));
    
    rtmin=fix(rt)-deltax+1;
        
    for i=1:N,      
        xtemp1=xtemp+rtmin(1,i)-1;
        ytemp1=ytemp+rtmin(2,i)-1;
        ztemp1=ztemp+rtmin(3,i)-1;
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=sub2ind_AS(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
        
        Im_S=zeros(S); 
        Im_S(temp_ind)=Im(indIm);
        
        Fxx=Im_S;
        Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2./(sigmax)^2;
        Fxtemp(:,i)=Fxx(:);
        
        Fyy=Im_S;
        Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2./(sigmay)^2;
        Fytemp(:,i)=Fyy(:);
        
        Fzz=Im_S;
        Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2./(sigmaz)^2;
        Fztemp(:,i)=Fzz(:);
        
        Itemp(:,i)=Im_S(:);           
    end

    Dg(1,:)=sum(G'.*Fxtemp);
    Dg(2,:)=sum(G'.*Fytemp);
    Dg(3,:)=sum(G'.*Fztemp);
    I_snake=sum(G'.*Itemp);

    vardelr2=rt*(diag(sum(AM,2))-AM);
    temp=Dg./(1+bpgrad.*BP);
    
    rt(:,notendp)=rt(:,notendp)+betta.*(temp(:,notendp)-2*alpha*vardelr2(:,notendp));
    lll(AM>0)=sum((rt(:,it)-rt(:,jt)).^2).^0.5;
    DgIav=sum(I_snake*diag(sum(lll))/sum(lll(:)));
    
    full([f1, N, sum(lll(:))/2, DgIav])
end

AM = LabelBranchesAM_AS(AM>0);
rt=rt';











% This version works with AM or AMlbl. Trees can be optimized seporately.
% Branch and end points can also be optimized
% Optimize_bps = 1,0 optimize branch points.
% Optimize_tps = 1,0 optimize terminal (start, end) points.
% Multiple_trees = 1,0 optimize labeled trees seporately.
% sig_ips, sig_bps, and sig_tps are the corresponding sigmas.

function [AMlbl rt I_snake]=Snake_All_test(Orig,AMlbl,r,Optimize_bps,Optimize_tps,Multiple_trees,pointsperum,Nstep,alpha,betta,sig_ips,sig_bps,sig_tps)

rt=r';
AM = spones(AMlbl+AMlbl');

if Multiple_trees==1
    % disconnect different label trees
    bp=find(sum(AM)>2);
    for i=1:length(bp)
        bp_labels=nonzeros(AMlbl(:,bp(i)));
        ubp_labels=unique(bp_labels);
        if length(ubp_labels)>1
            rt=[rt,rt(:,bp(i))*ones(1,length(ubp_labels)-1)];
            a=length(AM);
            AM(a+length(ubp_labels)-1,a+length(ubp_labels)-1)=0;
            for j=2:length(ubp_labels)
                temp=(AMlbl(:,bp(i))==ubp_labels(j));
                AM(temp,bp(i))=0;
                AM(bp(i),temp)=0;
                AM(temp,a+j-1)=1;
                AM(a+j-1,temp)=1;
            end
        end
    end
end

N=length(AM);
[it,jt]=find(AM);
lll=AM; lll(AM>0)=sum((rt(:,it)-rt(:,jt)).^2).^0.5;

format short g
full([N, sum(lll(:))/2])

%M=double(max(Orig(:)));
M=double(mean(Orig(Orig>0)));
sizeIm=size(Orig);

deltax=fix(2*sig_ips);
deltay=fix(2*sig_ips);
deltaz=fix(2*sig_ips);
S=[2*deltax,2*deltay,2*deltaz]+1;
[xtemp,ytemp,ztemp]=ind2sub(S,1:prod(S));

if Optimize_bps==1
    deltax_bp=fix(2*sig_bps);
    deltay_bp=fix(2*sig_bps);
    deltaz_bp=fix(2*sig_bps);
    S_bp=[2*deltax_bp,2*deltay_bp,2*deltaz_bp]+1;
    [xtemp_bp,ytemp_bp,ztemp_bp]=ind2sub(S_bp,1:prod(S_bp));
end

if Optimize_tps==1
    deltax_endp=fix(3*sig_tps);
    deltay_endp=fix(3*sig_tps);
    deltaz_endp=fix(3*sig_tps);
    S_endp=2.*[deltax_endp,deltay_endp,deltaz_endp]+1;
    [xtemp_endp,ytemp_endp,ztemp_endp]=ind2sub(S_endp,1:prod(S_endp));
end


% calculation of curvature, K
rt_dot=zeros(3,size(AM_connect,1));
for k=1:size(rt,2)
    temp=find(AM(k,:));
    if length(temp)==2
        rt_dot(:,k)=(rt(:,temp(2))-rt(:,temp(1)))./2;
    end
end
rt_2dot=rt*AM-2.*r_connect;
K=sum(cross(rt_dot,rt_2dot).^2,1).^0.5./sum(rt_dot.^2,1).^1.5;
K(sum(AM)~=2)=0;


for f1=1:Nstep
    % resegment
    [AM rt] = Divide_Segments_AS(AM,rt,pointsperum);
    [AM rt] = Merge_Segments_AS(AM,rt,pointsperum);
    [it,jt]=find(AM);
    N=length(AM);
    lll=AM; lll(AM>0)=sum((rt(:,it)-rt(:,jt)).^2).^0.5;  
    intermp=find(sum(AM)==2);
    endp=find(sum(AM)==1);
    bp=find(sum(AM)>2);
    
    Dg=zeros(3,N); I_snake=zeros(1,N);
    
    Fxtemp=zeros(prod(S),length(intermp)); Fytemp=zeros(prod(S),length(intermp)); Fztemp=zeros(prod(S),length(intermp));
    Itemp=zeros(prod(S),length(intermp)); 

    Gx=exp(-((rt(1,intermp)-round(rt(1,intermp)))'*ones(1,prod(S))+deltax+1-ones(length(intermp),1)*xtemp).^2./2./sig_ips^2);
    Gy=exp(-((rt(2,intermp)-round(rt(2,intermp)))'*ones(1,prod(S))+deltay+1-ones(length(intermp),1)*ytemp).^2./2./sig_ips^2);
    Gz=exp(-((rt(3,intermp)-round(rt(3,intermp)))'*ones(1,prod(S))+deltaz+1-ones(length(intermp),1)*ztemp).^2./2./sig_ips^2);
    G=Gx.*Gy.*Gz;
    G=G./(sum(G,2)*ones(1,prod(S)));
            
    for i=1:length(intermp),
        xtemp1=xtemp+round(rt(1,intermp(i)))-deltax-1;
        ytemp1=ytemp+round(rt(2,intermp(i)))-deltay-1;
        ztemp1=ztemp+round(rt(3,intermp(i)))-deltaz-1;
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=sub2ind_AS(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
        
        Im_S=zeros(S); 
        Im_S(temp_ind)=double(Orig(indIm))./M;
        
        Fxx=Im_S;
        Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2;
        Fxtemp(:,i)=Fxx(:);
        
        Fyy=Im_S;
        Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2;
        Fytemp(:,i)=Fyy(:);
        
        Fzz=Im_S;
        Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2;
        Fztemp(:,i)=Fzz(:);
        
        Itemp(:,i)=Im_S(:);     
    end
    
    Dg(1,intermp)=sum(G'.*Fxtemp);
    Dg(2,intermp)=sum(G'.*Fytemp);
    Dg(3,intermp)=sum(G'.*Fztemp);
    I_snake(intermp)=sum(G'.*Itemp);
    
    if Optimize_bps==1
        Fxtemp=zeros(prod(S_bp),length(bp)); Fytemp=zeros(prod(S_bp),length(bp)); Fztemp=zeros(prod(S_bp),length(bp));
        Itemp=zeros(prod(S_bp),length(bp));
        
        Gx=exp(-((rt(1,bp)-round(rt(1,bp)))'*ones(1,prod(S_bp))+deltax_bp+1-ones(length(bp),1)*xtemp_bp).^2./2./sig_bps^2);
        Gy=exp(-((rt(2,bp)-round(rt(2,bp)))'*ones(1,prod(S_bp))+deltay_bp+1-ones(length(bp),1)*ytemp_bp).^2./2./sig_bps^2);
        Gz=exp(-((rt(3,bp)-round(rt(3,bp)))'*ones(1,prod(S_bp))+deltaz_bp+1-ones(length(bp),1)*ztemp_bp).^2./2./sig_bps^2);
        G=Gx.*Gy.*Gz;
        G=G./(sum(G,2)*ones(1,prod(S_bp)));
        
        for i=1:length(bp),
            xtemp1=xtemp_bp+round(rt(1,bp(i)))-deltax_bp-1;
            ytemp1=ytemp_bp+round(rt(2,bp(i)))-deltay_bp-1;
            ztemp1=ztemp_bp+round(rt(3,bp(i)))-deltaz_bp-1;
            temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
            indIm=sub2ind_AS(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
            
            Im_S=zeros(S_bp);
            Im_S(temp_ind)=double(Orig(indIm))./M;
            
            Fxx=Im_S;
            Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2;
            Fxtemp(:,i)=Fxx(:);
            
            Fyy=Im_S;
            Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2;
            Fytemp(:,i)=Fyy(:);
            
            Fzz=Im_S;
            Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2;
            Fztemp(:,i)=Fzz(:);
            
            Itemp(:,i)=Im_S(:);
        end
        
        Dg(1,bp)=sum(G'.*Fxtemp);
        Dg(2,bp)=sum(G'.*Fytemp);
        Dg(3,bp)=sum(G'.*Fztemp);
        I_snake(bp)=sum(G'.*Itemp);
    end
    
    if Optimize_tps==1
        Fxtemp=zeros(prod(S_endp),length(endp)); Fytemp=zeros(prod(S_endp),length(endp)); Fztemp=zeros(prod(S_endp),length(endp));
        Itemp=zeros(prod(S_endp),length(endp));
        
        Gx=exp(-((rt(1,endp)-round(rt(1,endp)))'*ones(1,prod(S_endp))+deltax_endp+1-ones(length(endp),1)*xtemp_endp).^2./2./sig_tps^2);
        Gy=exp(-((rt(2,endp)-round(rt(2,endp)))'*ones(1,prod(S_endp))+deltay_endp+1-ones(length(endp),1)*ytemp_endp).^2./2./sig_tps^2);
        Gz=exp(-((rt(3,endp)-round(rt(3,endp)))'*ones(1,prod(S_endp))+deltaz_endp+1-ones(length(endp),1)*ztemp_endp).^2./2./sig_tps^2);
        G=Gx.*Gy.*Gz;
        G=G./(sum(G,2)*ones(1,prod(S_endp)));
        G=G.*cos(((((rt(1,endp)-round(rt(1,endp)))'*ones(1,prod(S_endp))+deltax_endp+1-ones(length(endp),1)*xtemp_endp)./sig_tps).^2+(((rt(2,endp)-round(rt(2,endp)))'*ones(1,prod(S_endp))+deltay_endp+1-ones(length(endp),1)*ytemp_endp)./sig_tps).^2+(((rt(3,endp)-round(rt(3,endp)))'*ones(1,prod(S_endp))+deltaz_endp+1-ones(length(endp),1)*ztemp_endp)./sig_tps).^2).^0.5);
        
        for i=1:length(endp),
            xtemp1=xtemp_endp+round(rt(1,endp(i)))-deltax_endp-1;
            ytemp1=ytemp_endp+round(rt(2,endp(i)))-deltay_endp-1;
            ztemp1=ztemp_endp+round(rt(3,endp(i)))-deltaz_endp-1;
            temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
            indIm=sub2ind_AS(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
            
            Im_S=zeros(S_endp);
            Im_S(temp_ind)=double(Orig(indIm))./M;
            
            Fxx=Im_S;
            Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2;
            Fxtemp(:,i)=Fxx(:);
            
            Fyy=Im_S;
            Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2;
            Fytemp(:,i)=Fyy(:);
            
            Fzz=Im_S;
            Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2;
            Fztemp(:,i)=Fzz(:);
            
            Itemp(:,i)=Im_S(:);
        end
        
        Dg(1,endp)=sum(G'.*Fxtemp);
        Dg(2,endp)=sum(G'.*Fytemp);
        Dg(3,endp)=sum(G'.*Fztemp);
        I_snake(endp)=sum(G'.*Itemp);
    end
    
    vardelr2=rt*(diag(sum(AM,2))-AM);
    
    rt(:,intermp)=rt(:,intermp)+betta.*(Dg(:,intermp)-2*alpha*vardelr2(:,intermp));
    if Optimize_bps==1
        rt(:,bp)=rt(:,bp)+betta.*(100.*Dg(:,bp)-2*alpha*vardelr2(:,bp));
    end
    if Optimize_tps==1
        rt(:,endp)=rt(:,endp)+betta.*(100.*Dg(:,endp)-2*alpha*vardelr2(:,endp));
    end
    
    lll(AM>0)=sum((rt(:,it)-rt(:,jt)).^2).^0.5;
    DgIav=sum(I_snake*diag(sum(lll))/sum(lll(:)));
    
    full([f1, N, sum(lll(:))/2, DgIav])
end

if Multiple_trees==0
    AMlbl = LabelBranchesAM_AS(AM>0);
elseif Multiple_trees==1
    AMlbl = LabelTreesAM_AS(AM);
end

rt=rt';











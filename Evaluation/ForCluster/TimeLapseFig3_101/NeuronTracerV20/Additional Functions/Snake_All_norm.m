% This version works with AM or AMlbl. Trees can be optimized separately.
% Branch and end points can also be optimized
% Optimize_bps = 1,0 optimize branch points.
% Optimize_tps = 1,0 optimize terminal (start, end) points.
% Multiple_trees = 1,0 optimize labeled trees separately.
% sig_ips, sig_bps, and sig_tps are the corresponding sigmas.
% this version of the code is normalized for ppm as in the paper

function [AMlbl r R I_snake]=Snake_All_norm(Orig,AMlbl,r,R,Optimize_bps,Optimize_tps,Multiple_trees,pointsperum,Nstep,alpha,betta,sig_ips,sig_bps,sig_tps,output)

if output==1
    disp('F1 trace optimization started.')
    format short g
    display(['      Iteration   ',   '# vertices   ',   'Trace length   ',    '<Intensity>'])
end

AMlbl=max(AMlbl,AMlbl');
rem_ind=(sum(AMlbl,1)==0);
AMlbl(:,rem_ind)=[];
AMlbl(rem_ind,:)=[];
r(rem_ind,:)=[];
R(rem_ind)=[];

AM = spones(AMlbl+AMlbl');

if Multiple_trees==1
    % disconnect different label trees
    bp=find(sum(AM)>2);
    for i=1:length(bp)
        bp_labels=nonzeros(AMlbl(:,bp(i)));
        ubp_labels=unique(bp_labels);
        if length(ubp_labels)>1
            r=[r;ones(length(ubp_labels)-1,1)*r(bp(i),:)];
            R=[R;ones(length(ubp_labels)-1,1).*R(bp(i))];
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

[it,jt]=find(AM);
lll=AM; lll(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^0.5;
L0=sum(lll(:))/2;

M=double(max(max(max(Orig))));
sizeIm=size(Orig);

deltax=fix(2*sig_ips);
deltay=fix(2*sig_ips);
deltaz=fix(2*sig_ips);
S=[2*deltax,2*deltay,2*deltaz]+1;
[xtemp,ytemp,ztemp]=ind2sub(S,1:prod(S));

deltax_bp=fix(2*sig_bps);
deltay_bp=fix(2*sig_bps);
deltaz_bp=fix(2*sig_bps);
S_bp=[2*deltax_bp,2*deltay_bp,2*deltaz_bp]+1;
[xtemp_bp,ytemp_bp,ztemp_bp]=ind2sub(S_bp,1:prod(S_bp));

deltax_endp=fix(3*sig_tps);
deltay_endp=fix(3*sig_tps);
deltaz_endp=fix(3*sig_tps);
S_endp=2.*[deltax_endp,deltay_endp,deltaz_endp]+1;
[xtemp_endp,ytemp_endp,ztemp_endp]=ind2sub(S_endp,1:prod(S_endp));

for f1=1:Nstep
    % resegment
    [AM r R] = Divide_Segments(AM,r,R,pointsperum);
    [AM r R] = Merge_Segments(AM,r,R,pointsperum);
    [it,jt]=find(AM);
    N=length(AM);
    lll=AM; lll(AM>0)=sum((r(it,:)-r(jt,:)).^2,2).^0.5; 
    L=sum(lll(:))/2;
    if L/L0>3
        error('Trace is unstable. Reduce Trace Stiffness (alpha) or Optimization Step Size (beta).')
    end
    
    intermp=find(sum(AM)==2);
    endp=find(sum(AM)==1);
    bp=find(sum(AM)>2);
    
    Dg=zeros(N,3); I_snake=zeros(N,1);
    
    Fxtemp=zeros(length(intermp),prod(S)); Fytemp=zeros(length(intermp),prod(S)); Fztemp=zeros(length(intermp),prod(S));
    Itemp=zeros(length(intermp),prod(S)); 

    Gx=exp(-((r(intermp,1)-round(r(intermp,1)))*ones(1,prod(S))+deltax+1-ones(length(intermp),1)*xtemp).^2./2./sig_ips^2);
    Gy=exp(-((r(intermp,2)-round(r(intermp,2)))*ones(1,prod(S))+deltay+1-ones(length(intermp),1)*ytemp).^2./2./sig_ips^2);
    Gz=exp(-((r(intermp,3)-round(r(intermp,3)))*ones(1,prod(S))+deltaz+1-ones(length(intermp),1)*ztemp).^2./2./sig_ips^2);
    G=Gx.*Gy.*Gz;
    G=G./(sum(G,2)*ones(1,prod(S)));
            
    for i=1:length(intermp),
        xtemp1=xtemp+round(r(intermp(i),1))-deltax-1;
        ytemp1=ytemp+round(r(intermp(i),2))-deltay-1;
        ztemp1=ztemp+round(r(intermp(i),3))-deltaz-1;
        temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
        indIm=sub2ind_ASfast(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
        
        Im_S=zeros(S); 
        Im_S(temp_ind)=double(Orig(indIm))./M;
        
        Fxx=Im_S;
        Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2;
        Fxtemp(i,:)=Fxx(:);
        
        Fyy=Im_S;
        Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2;
        Fytemp(i,:)=Fyy(:);
        
        Fzz=Im_S;
        Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2;
        Fztemp(i,:)=Fzz(:);
        
        Itemp(i,:)=Im_S(:);     
    end
    
    Dg(intermp,1)=sum(G.*Fxtemp,2);
    Dg(intermp,2)=sum(G.*Fytemp,2);
    Dg(intermp,3)=sum(G.*Fztemp,2);
    I_snake(intermp)=sum(G.*Itemp,2);
    
    if Optimize_bps==1
        Fxtemp=zeros(length(bp),prod(S_bp)); Fytemp=zeros(length(bp),prod(S_bp)); Fztemp=zeros(length(bp),prod(S_bp));
        Itemp=zeros(length(bp),prod(S_bp));
        
        Gx=exp(-((r(bp,1)-round(r(bp,1)))*ones(1,prod(S_bp))+deltax_bp+1-ones(length(bp),1)*xtemp_bp).^2./2./sig_bps^2);
        Gy=exp(-((r(bp,2)-round(r(bp,2)))*ones(1,prod(S_bp))+deltay_bp+1-ones(length(bp),1)*ytemp_bp).^2./2./sig_bps^2);
        Gz=exp(-((r(bp,3)-round(r(bp,3)))*ones(1,prod(S_bp))+deltaz_bp+1-ones(length(bp),1)*ztemp_bp).^2./2./sig_bps^2);
        G=Gx.*Gy.*Gz;
        G=G./(sum(G,2)*ones(1,prod(S_bp)));
        
        for i=1:length(bp),
            xtemp1=xtemp_bp+round(r(bp(i),1))-deltax_bp-1;
            ytemp1=ytemp_bp+round(r(bp(i),2))-deltay_bp-1;
            ztemp1=ztemp_bp+round(r(bp(i),3))-deltaz_bp-1;
            temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
            indIm=sub2ind_ASfast(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
            
            Im_S=zeros(S_bp);
            Im_S(temp_ind)=double(Orig(indIm))./M;
            
            Fxx=Im_S;
            Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2;
            Fxtemp(i,:)=Fxx(:);
            
            Fyy=Im_S;
            Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2;
            Fytemp(i,:)=Fyy(:);
            
            Fzz=Im_S;
            Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2;
            Fztemp(i,:)=Fzz(:);
            
            Itemp(i,:)=Im_S(:);
        end
        
        Dg(bp,1)=sum(G.*Fxtemp,2);
        Dg(bp,2)=sum(G.*Fytemp,2);
        Dg(bp,3)=sum(G.*Fztemp,2);
        I_snake(bp)=sum(G.*Itemp,2);
    end
    
    if Optimize_tps==1
        Fxtemp=zeros(length(endp),prod(S_endp)); Fytemp=zeros(length(endp),prod(S_endp)); Fztemp=zeros(length(endp),prod(S_endp));
        Itemp=zeros(length(endp),prod(S_endp));
        
        Gx=exp(-((r(endp,1)-round(r(endp,1)))*ones(1,prod(S_endp))+deltax_endp+1-ones(length(endp),1)*xtemp_endp).^2./2./sig_tps^2);
        Gy=exp(-((r(endp,2)-round(r(endp,2)))*ones(1,prod(S_endp))+deltay_endp+1-ones(length(endp),1)*ytemp_endp).^2./2./sig_tps^2);
        Gz=exp(-((r(endp,3)-round(r(endp,3)))*ones(1,prod(S_endp))+deltaz_endp+1-ones(length(endp),1)*ztemp_endp).^2./2./sig_tps^2);
        G=Gx.*Gy.*Gz;
        G=G./(sum(G,2)*ones(1,prod(S_endp)));
        G=G.*cos(((((r(endp,1)-round(r(endp,1)))*ones(1,prod(S_endp))+deltax_endp+1-ones(length(endp),1)*xtemp_endp)./sig_tps).^2+...
            (((r(endp,2)-round(r(endp,2)))*ones(1,prod(S_endp))+deltay_endp+1-ones(length(endp),1)*ytemp_endp)./sig_tps).^2+...
            (((r(endp,3)-round(r(endp,3)))*ones(1,prod(S_endp))+deltaz_endp+1-ones(length(endp),1)*ztemp_endp)./sig_tps).^2).^0.5);
        
        for i=1:length(endp),
            xtemp1=xtemp_endp+round(r(endp(i),1))-deltax_endp-1;
            ytemp1=ytemp_endp+round(r(endp(i),2))-deltay_endp-1;
            ztemp1=ztemp_endp+round(r(endp(i),3))-deltaz_endp-1;
            temp_ind=(xtemp1>=1 & xtemp1<=sizeIm(1) & ytemp1>=1 & ytemp1<=sizeIm(2) & ztemp1>=1 & ztemp1<=sizeIm(3));
            indIm=sub2ind_ASfast(sizeIm,xtemp1(temp_ind),ytemp1(temp_ind),ztemp1(temp_ind));
            
            Im_S=zeros(S_endp);
            Im_S(temp_ind)=double(Orig(indIm))./M;
            
            Fxx=Im_S;
            Fxx(2:end-1,:,:)=(Fxx(3:end,:,:)-Fxx(1:end-2,:,:))./2;
            Fxtemp(i,:)=Fxx(:);
            
            Fyy=Im_S;
            Fyy(:,2:end-1,:)=(Fyy(:,3:end,:)-Fyy(:,1:end-2,:))./2;
            Fytemp(i,:)=Fyy(:);
            
            Fzz=Im_S;
            Fzz(:,:,2:end-1)=(Fzz(:,:,3:end)-Fzz(:,:,1:end-2))./2;
            Fztemp(i,:)=Fzz(:);
            
            Itemp(i,:)=Im_S(:);
        end
        
        Dg(endp,1)=sum(G.*Fxtemp,2);
        Dg(endp,2)=sum(G.*Fytemp,2);
        Dg(endp,3)=sum(G.*Fztemp,2);
        I_snake(endp)=sum(G.*Itemp,2);
    end
    
    if output==1
        DgIav=sum(diag(sum(lll))/sum(lll(:))*I_snake);
        disp(full([f1, N, L, DgIav]))
    end
    
    vardelr2=(diag(sum(AM,2))-AM)*r;
    
    r(intermp,:)=r(intermp,:)+betta.*(Dg(intermp,:)./pointsperum-2*alpha.*pointsperum*vardelr2(intermp,:));
    if Optimize_bps==1
        r(bp,:)=r(bp,:)+betta.*((sig_bps/sig_ips)^3.*Dg(bp,:)./pointsperum-2*alpha.*pointsperum*vardelr2(bp,:));
    end
    if Optimize_tps==1
        r(endp,:)=r(endp,:)+betta.*((sig_tps/sig_ips)^3.*Dg(endp,:)./pointsperum-2*alpha.*pointsperum*vardelr2(endp,:));
    end
end

if Multiple_trees==0
    AMlbl = LabelBranchesAM(AM>0);
elseif Multiple_trees==1
    AMlbl = LabelTreesAM(AM);
end

if output==1
    disp('F1 trace optimization is complete.')
end









% This function returns a=true if two parallelepipeds overlap

function a=Overlap(Verts1,Verts2)

a=true;
faceind=[1,2,4,3;5,6,8,7;3,4,8,7;1,2,6,5;1,3,7,5;2,4,8,6];

i=1;
while a && i<=6
    normal=cross(Verts1(faceind(i,2),:)-Verts1(faceind(i,1),:),Verts1(faceind(i,4),:)-Verts1(faceind(i,1),:));
    other_verts=1:8;
    other_verts(faceind(i,:))=[];
    side1=(Verts1(other_verts(1),:)-Verts1(faceind(i,1),:))*normal';
    side2=(Verts2-ones(8,1)*Verts1(faceind(i,1),:))*normal';
    if nnz(side1.*side2>0)==0
        a=false;
    end
    i=i+1;
end

i=1;
while a && i<=6
    normal=cross(Verts2(faceind(i,2),:)-Verts2(faceind(i,1),:),Verts2(faceind(i,4),:)-Verts2(faceind(i,1),:));
    other_verts=1:8;
    other_verts(faceind(i,:))=[];
    side1=(Verts2(other_verts(1),:)-Verts2(faceind(i,1),:))*normal';
    side2=(Verts1-ones(8,1)*Verts2(faceind(i,1),:))*normal';
    if nnz(side1.*side2>0)==0
        a=false;
    end
    i=i+1;
end
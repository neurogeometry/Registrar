function createXML(TilePositions,TileSize_org,SaveFolder)


docNode = com.mathworks.xml.XMLUtils.createDocument('SpimData')
toc = docNode.getDocumentElement;
toc.setAttribute('version','0.2');

product = docNode.createElement('BasePath');
product.setAttribute('type','relative');
product.appendChild(docNode.createTextNode('.'));
toc.appendChild(product)

product = docNode.createElement('SequenceDescription');
toc.appendChild(product)



    curr_node = docNode.createElement('ImageLoader');
    curr_node.setAttribute('format','bdv.hdf5');
    product.appendChild(curr_node);
    
    
    curr_node1 = docNode.createElement('hdf5');
    curr_node1.setAttribute('type','relative');
    curr_node1.appendChild(docNode.createTextNode('export1.h5'));
    curr_node.appendChild(curr_node1);
    
    curr_node1 = docNode.createElement('ViewSetups');
    product.appendChild(curr_node1);


for idx = 0:size(TilePositions,1)-1
    curr_node = docNode.createElement('ViewSetup');
    curr_node1.appendChild(curr_node);
    
        curr_node2 = docNode.createElement('id');
        curr_node2.appendChild(docNode.createTextNode(num2str(idx)));
        curr_node.appendChild(curr_node2);
        
        curr_node2 = docNode.createElement('name');
        curr_node2.appendChild(docNode.createTextNode('channel 1'));
        curr_node.appendChild(curr_node2);
        
        curr_node2 = docNode.createElement('size');
        curr_node2.appendChild(docNode.createTextNode([num2str(TileSize_org(1)),' ',num2str(TileSize_org(2)),' ',num2str(TileSize_org(3))]));
        curr_node.appendChild(curr_node2);
        
        curr_node2 = docNode.createElement('voxelSize');
        curr_node.appendChild(curr_node2);
        
            curr_node3 = docNode.createElement('unit');
            curr_node3.appendChild(docNode.createTextNode('pixel'));
            curr_node2.appendChild(curr_node3);
            
            curr_node3 = docNode.createElement('size');
            curr_node3.appendChild(docNode.createTextNode('1.0 1.0 1.0'));
            curr_node2.appendChild(curr_node3);
            
        curr_node2 = docNode.createElement('attributes');
        curr_node.appendChild(curr_node2);
        
            curr_node3 = docNode.createElement('channel');
            curr_node3.appendChild(docNode.createTextNode('1'));
            curr_node2.appendChild(curr_node3);
    
    curr_node = docNode.createElement('Attributes');
    curr_node.setAttribute('name','channel');
    curr_node1.appendChild(curr_node);
    
        curr_node2 = docNode.createElement('Channel');
        curr_node.appendChild(curr_node2);
        
            curr_node3 = docNode.createElement('id');
            curr_node3.appendChild(docNode.createTextNode('1'));
            curr_node2.appendChild(curr_node3);
            
            curr_node3 = docNode.createElement('name');
            curr_node3.appendChild(docNode.createTextNode('1'));
            curr_node2.appendChild(curr_node3);
end

    curr_node1 = docNode.createElement('Timepoints');
    curr_node1.setAttribute('type','range');
    product.appendChild(curr_node1);
    
    curr_node = docNode.createElement('first');
    curr_node.appendChild(docNode.createTextNode('0'));
    curr_node1.appendChild(curr_node);
    
    curr_node = docNode.createElement('last');
    curr_node.appendChild(docNode.createTextNode('0'));
    curr_node1.appendChild(curr_node);
    
    
product = docNode.createElement('ViewRegistrations');
toc.appendChild(product)
for idx = 0:size(TilePositions,1)-1
    
    
    x = TilePositions(idx+1,1);
    y = TilePositions(idx+1,2);
    z = TilePositions(idx+1,3);
    
    curr_node1 = docNode.createElement('ViewRegistration');
    curr_node1.setAttribute('setup',num2str(idx));
    curr_node1.setAttribute('timepoint',num2str(0));
    product.appendChild(curr_node1);
    
        curr_node2 = docNode.createElement('ViewTransform');
        curr_node2.setAttribute('type','affine');
        curr_node1.appendChild(curr_node2);
        
            curr_node3 = docNode.createElement('affine');
            curr_node3.appendChild(docNode.createTextNode(['1.0 0.0 0.0 ',num2str(x),' 0.0 1.0 0.0 ',num2str(y),' 0.0 0.0 1.0 ',num2str(z)]));
            curr_node2.appendChild(curr_node3);
end
xmlFileName = [SaveFolder,'HDf5\export1.xml'];
xmlwrite(xmlFileName,docNode);
type(xmlFileName);
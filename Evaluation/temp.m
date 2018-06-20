n_points = 20; 
  source = 10 * rand(n_points, 2); 
  target = source + rand(n_points, 2); 
  target_indices = nearestneighborlinker(source, target); 
  colors = hsv(n_points); 
  figure 
  hold on 
  for i = 1 :n_points 
     plot(source(i,1), source(i,2), 'o', 'Color', colors(i,:)) 
     plot(target(target_indices(i),1), target(target_indices(i),2), 's', ... 
        'Color', colors(i,:)) 
     plot( [ source(i,1) target(target_indices(i),1) ] , ... 
        [ source(i,2) target(target_indices(i),2) ], ... 
         'Color', colors(i,:)) 
  end 
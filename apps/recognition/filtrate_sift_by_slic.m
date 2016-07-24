function feature=filtrate_sift_by_slic(hight,s,features)
%filtrate sift descriptor by slic edge. 
%%wight is the wight of image
%hight is the hight of image
%s is the matrix include the edge of slic
%features is the sift features

for i = 1:size(s,1)
    temp(1,i) = mod(s(i,1) , hight);
    temp(2,i) = fix(s(i,1) / hight) + 1;
end

% Buld a kd-tree
kdtree = vl_kdtreebuild( temp ) ;

Q = [features.frame(1,1);features.frame(2,1)];
[index, dist] = vl_kdtreequery (kdtree, temp, Q, 'MaxNumComparisons', 1);
if dist <= features.frame(3,1);
    features.frame(:,1) = 0;
    features.descr(:,1) = 0;
    features.contrast(:,1) = 0;
end

for i = 2:size(features.frame,2)    
    
    a = 8 / 3 * abs(temp(1,index) - features.frame(1,i));
    a = a * a;
    b = 8 / 3 * abs(temp(2,index) - features.frame(2,i));
    b = b * b;
    c = features.frame(3,i) * features.frame(3,i);
    if a + b <= c
        features.frame(:,i) = 0;
        features.descr(:,i) = 0;
        features.contrast(:,i) = 0;
        continue
        
    else
        Q = [features.frame(1,i);features.frame(2,i)];
        [index, dist] = vl_kdtreequery (kdtree, temp, Q, 'MaxNumComparisons', 1);
        if dist <= features.frame(3,i);
            features.frame(:,i) = 0;
            features.descr(:,i) = 0;
            features.contrast(:,i) = 0;
        end
    end
    
end

features.frame(:,all(features.frame==0,1))=[] ;
features.descr(:,all(features.descr==0,1))=[] ;
features.contrast(:,all(features.contrast==0,1))=[] ;

feature = features;
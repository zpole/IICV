function feature=filtrate_sift_by_slic(hight,s,features, opts)
%filtrate sift descriptor by slic edge. 
%%wight is the wight of image
%hight is the hight of image
%s is the matrix include the edge of slic
%features is the sift features for dsift,features must include three matrix
%frame, descr and contrast(look vl_dsift.m); for densecolor, features
%include two matrix dcolor and info dcolor.
%opts have two options now, 'dsift' and 'dcolor'.

for i = 1:size(s,1)
    temp(2,i) = mod(s(i,1) , hight);
    temp(1,i) = fix(s(i,1) / hight) + 1;
end

% Buld a kd-tree
kdtree = vl_kdtreebuild( temp ) ;

switch opts
    case 'dsift'
        sel = vl_colsubset(1:size(features.descr,2), 5120) ;
        features.frame = features.frame(:,sel) ;
        features.descr = features.descr(:,sel) ;
        features.contrast = features.contrast(:,sel) ;
        fprintf('%d,%d\n',size(s,1),size(features.frame,2));
        %Q = [features.frame(1,1);features.frame(2,1)];
        %[index, dist] = vl_kdtreequery (kdtree, temp, Q);
        %if dist <= features.frame(3,1);
        %    features.frame(:,1) = 0;
        %    features.descr(:,1) = 0;
        %    features.contrast(:,1) = 0;
        %end

        for i = 1:size(features.frame,2)    
    
        %    a = 4 / 3 * abs(temp(1,index) - features.frame(1,i));
        %    a = a * a;
        %    b = 4 / 3 * abs(temp(2,index) - features.frame(2,i));
        %    b = b * b;
        %    c = features.frame(3,i) * features.frame(3,i);
        %    if a + b <= c
        %        features.frame(:,i) = 0;
        %        features.descr(:,i) = 0;
        %        features.contrast(:,i) = 0;
        %        continue
        
       %     else
                Q = [features.frame(1,i);features.frame(2,i)];
                [index, dist] = vl_kdtreequery (kdtree, temp, Q);
                if dist > features.frame(3,i);
                    features.frame(:,i) = -1;
                    features.descr(:,i) = -1;
                    features.contrast(:,i) = -1;
                end
        %    end
    
        end

        features.frame(:,all(features.frame==-1,1))=[] ;
        features.descr(:,all(features.descr==-1,1))=[] ;
        features.contrast(:,all(features.contrast==-1,1))=[] ;
        
    case 'dcolor'
        for i = 1:size(features.frame,2)    
                Q = [features.frame(1,i);features.frame(2,i)];
                [index, dist] = vl_kdtreequery (kdtree, temp, Q);
                if dist > features.frame(3,i);
                    features.frame(:,i) = -1;
                    features.descr(:,i) = -1;
                end
        end

        features.frame(:,all(features.frame==-1,1))=[] ;
        features.descr(:,all(features.descr==-1,1))=[] ;
end

feature = features;
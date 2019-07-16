% Rough estimates for principal point and focal length:
% fg = [1009.661057548567   1009.706970843858 ]';
% cg = [ 1031.688786545889   1288.395050759215]';
% kg = [-0.054285891395760  -0.026812583950935 0 0]';

% % wuyuan fisheye test
% fg = [320   280 ]';
% cg = [ 640   360]';
% kg = [0.2  -0.1 0.04 -0.008]';

%opencv init
% fg = [350.438926029233   350.2174211053091 ]';
% cg = [643.9272197363401   359.684099539598]';
% kg = [0.121311, -0.0902425, 0.10761, -0.0461907]';

% fg = [320.88785   320.48022 ]';
% cg = [656.05855   379.36670]';
% kg = [0.06336   0.06080   -0.03398   0.00450]';

%%%最近使用的参数
fg = [377.52445   363.102738465345680  ]';
cg = [637.656119618039840  367.716794364937930 ]';
kg = [-0.052199196426168  0.020099287047905  -0.018613380783400  0.004923019763423]';

grid_success = 0;


while (~grid_success)
    
    if  1%manual_squares,
        fprintf(1,['Number of squares along the X direction ([]=' num2str(n_sq_x_default) ') = ']);
        n_sq_x = input(''); %6

        if isempty(n_sq_x), n_sq_x = n_sq_x_default; end;
        fprintf(1,['Number of squares along the Y direction ([]=' num2str(n_sq_y_default) ') = ']);
        n_sq_y = input(''); %6

        if isempty(n_sq_y), n_sq_y = n_sq_y_default; end;
        Np = (n_sq_x+1)*(n_sq_y+1);
        eval(['x_' num2str(kk) ' =  zeros(2,' num2str(Np) ') ;']);
        
        grid_success = 1;

        figure(2); clf;
        image(uint8(I));
        axis image;
        colormap(map);
        set(2,'color',[1 1 1]);
        x=zeros(2,Np);y = zeros(2,Np);
        figure(2); hold on;
        v0=0;
        [X,Y,cropimg,RECT] = imcrop(uint8(I));
        corners = findCorners(cropimg,0.05,1);
        chessboards = chessboardsFromCorners(corners);
        
        x_fc = zeros(4,1);
        y_fc = zeros(4,1);
    
        for i=1:length(chessboards)
            % extract chessboard
            cb = chessboards{i};
            A=zeros(size(cb,1),size(cb,2),2);
            for j=1:size(cb,1)
                p = corners.p(cb(j,:),:)+1;
                %plot(RECT(1)+p(:,1),RECT(2)+v0+p(:,2),'o','Color','r','LineWidth',2);
                A(j,:,:)=[RECT(1)+p(:,1),RECT(2)+v0+p(:,2)];
            end
        end
        ns=size(A);
        if(ns(1)~=n_sq_y_default+1)
            A=permute(A,[2 1 3]);
            A=flipud(A);
        end
        for j=1:size(A,1)
            plot(A(j,:,1),A(j,:,2),'o','Color','r','LineWidth',2);
        end
        [hang,lie,ch]=size(A);
        if(hang*lie~=Np)
            error('detect chessboard size error!');
        end
        for j=1:Np
            if(mod(j,lie)==0)
                hi=floor(j/lie);
                hj=lie;
            else
                hi=floor(j/lie)+1;
                hj=mod(j,lie);
            end
            text(A(hi,hj,1),A(hi,hj,2)+20,num2str(j),'horiz','center','color','r');
        end
        hold off;
        figure(2); 
        image(I);
        colormap(map);
        hold on;

        hx = plot(zeros(1,Np),zeros(1,Np),'r+');
        hcp = plot(1,1,'co');

        hold off;
        
        for np = 1:Np
            set(hcp,'Xdata',x(1,np)+1,'Ydata',x(2,np)+1);

            % [xi,yi,b] = ginput4(1);
            if(mod(np,lie)==0)
                hi=floor(np/lie);
                hj=lie;
            else
                hi=floor(np/lie)+1;
                hj=mod(np,lie);
            end
            b=1;
            xi=A(hi,hj,1);
            yi=A(hi,hj,2);
            if b==1,
                xxi = [xi;yi];
                if q_converge,
                    [xxi] = cornerfinder(xxi,I,winty,wintx);
                end;
                x(1,np) = xxi(1) - 1;
                x(2,np) = xxi(2) - 1;
                set(hx,'Xdata',x(1,:)+1,'Ydata',x(2,:)+1);
                x_fc =[x(1,end-n_sq_x),x(1,1) ,x(1,1+n_sq_x),x(1,end)]';
                y_fc =[x(2,end-n_sq_x),x(2,1) ,x(2,1+n_sq_x),x(2,end)]';    
            end;
        end;
    end
    
    
    
    % Sort the corners:
    x_mean = mean(x_fc);
    y_mean = mean(y_fc);
    x_v = x_fc - x_mean;
    y_v = y_fc - y_mean;

    theta = atan2(-y_v,x_v);
    [junk,ind] = sort(theta);

    [junk,ind] = sort(mod(theta-theta(1),2*pi));

    ind = ind([4 3 2 1]); %-> New: the Z axis is pointing uppward

    x = x_fc(ind);
    y = y_fc(ind);

    x1= x(1); x2 = x(2); x3 = x(3); x4 = x(4);
    y1= y(1); y2 = y(2); y3 = y(3); y4 = y(4);

    % Find center:
    p_center = cross(cross([x1;y1;1],[x3;y3;1]),cross([x2;y2;1],[x4;y4;1]));
    x5 = p_center(1)/p_center(3);
    y5 = p_center(2)/p_center(3);

    % center on the X axis:
    x6 = (x3 + x4)/2;
    y6 = (y3 + y4)/2;

    % center on the Y axis:
    x7 = (x1 + x4)/2;
    y7 = (y1 + y4)/2;

    % Direction of displacement for the X axis:
    vX = [x6-x5;y6-y5];
    vX = vX / norm(vX);

    % Direction of displacement for the X axis:
    vY = [x7-x5;y7-y5];
    vY = vY / norm(vY);

    % Direction of diagonal:
    vO = [x4 - x5; y4 - y5];
    vO = vO / norm(vO);

    delta = 30;
    
    if ~grid_success
        fprintf(1,'Invalid grid. Try again.\n');
    end;
end;

n_sq_x_default = n_sq_x;
n_sq_y_default = n_sq_y;

x_n = (x - 1 - cg(1))/fg(1);
y_n = (y - 1 - cg(2))/fg(2);

[x_pn] = comp_fisheye_distortion([x_n' ; y_n'],kg);

% Compute the inside points through computation of the planar homography (collineation)
a00 = [x_pn(1,1);x_pn(2,1);1];
a10 = [x_pn(1,2);x_pn(2,2);1];
a11 = [x_pn(1,3);x_pn(2,3);1];
a01 = [x_pn(1,4);x_pn(2,4);1];

% Compute the planar collineation: (return the normalization matrix as well)
[Homo,Hnorm,inv_Hnorm] = compute_homography([a00 a10 a11 a01],[0 1 1 0;0 0 1 1;1 1 1 1]);

% Build the grid using the planar collineation:
x_l = ((0:n_sq_x)'*ones(1,n_sq_y+1))/n_sq_x;
y_l = (ones(n_sq_x+1,1)*(0:n_sq_y))/n_sq_y;
pts = [x_l(:) y_l(:) ones((n_sq_x+1)*(n_sq_y+1),1)]';

XXpn = Homo*pts;
XXpn = XXpn(1:2,:) ./ (ones(2,1)*XXpn(3,:));

XX = apply_fisheye_distortion(XXpn,kg);

XX(1,:) = fg(1)*XX(1,:) + cg(1) + 1;
XX(2,:) = fg(2)*XX(2,:) + cg(2) + 1;

% Complete size of the rectangle
W = n_sq_x*dX;
L = n_sq_y*dY;

Np = (n_sq_x+1)*(n_sq_y+1);
disp('Corner extraction...');
grid_pts = cornerfinder(XX,I,winty,wintx); %%% Finds the exact corners at every points!
%grid_pts = XX; %%% Finds the exact corners at every points!

grid_pts = grid_pts - 1; % subtract 1 to bring the origin to (0,0) instead of (1,1) in matlab (not necessary in C)

ind_corners = [1 n_sq_x+1 (n_sq_x+1)*n_sq_y+1 (n_sq_x+1)*(n_sq_y+1)]; % index of the 4 corners
ind_orig = (n_sq_x+1)*n_sq_y + 1;
xorig = grid_pts(1,ind_orig);
yorig = grid_pts(2,ind_orig);
dxpos = mean([grid_pts(:,ind_orig) grid_pts(:,ind_orig+1)]');
dypos = mean([grid_pts(:,ind_orig) grid_pts(:,ind_orig-n_sq_x-1)]');

x_box_kk = [grid_pts(1,:)-(wintx+.5);grid_pts(1,:)+(wintx+.5);grid_pts(1,:)+(wintx+.5);grid_pts(1,:)-(wintx+.5);grid_pts(1,:)-(wintx+.5)];
y_box_kk = [grid_pts(2,:)-(winty+.5);grid_pts(2,:)-(winty+.5);grid_pts(2,:)+(winty+.5);grid_pts(2,:)+(winty+.5);grid_pts(2,:)-(winty+.5)];


Xi = reshape(([0:n_sq_x]*dX)'*ones(1,n_sq_y+1),Np,1)';
Yi = reshape(ones(n_sq_x+1,1)*[n_sq_y:-1:0]*dY,Np,1)';
Zi = zeros(1,Np);

Xgrid = [Xi;Yi;Zi];

% All the point coordinates (on the image, and in 3D) - for global optimization:
x = grid_pts;
X = Xgrid;

%Saves all the data into variables:
eval(['wintx_' num2str(kk) ' = wintx;']);
eval(['winty_' num2str(kk) ' = winty;']);
eval(['dX_' num2str(kk) ' = dX;']);
eval(['dY_' num2str(kk) ' = dY;']);  
eval(['x_' num2str(kk) ' = x;']);%x_1=x;
eval(['X_' num2str(kk) ' = X;']);
eval(['n_sq_x_' num2str(kk) ' = n_sq_x;']);
eval(['n_sq_y_' num2str(kk) ' = n_sq_y;']);

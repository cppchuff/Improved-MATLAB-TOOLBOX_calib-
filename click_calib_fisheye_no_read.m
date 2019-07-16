%if exist('images_read');
%   active_images = active_images & images_read;
%end;

var2fix = 'dX_default';

fixvariable;

var2fix = 'dY_default';

fixvariable;

var2fix = 'map';

fixvariable;


if ~exist('n_ima'),
    data_calib_no_read;
end;

check_active_images;

% Step used to clean the memory if a previous atttempt has been made to read the entire set of images into memory:
for kk = 1:n_ima,
    if (exist(['I_' num2str(kk)])==1),
        clear(['I_' num2str(kk)]);
    end;
end;

fprintf(1,'\nManual re-extraction of the grid corners on the images\n');
fprintf(1,'Do you want to try to automatically find the closest corner? - only works with ckecker board corners  ([]=yes, other = no)\n');
q_converge = input('');

if isempty(q_converge),
    q_converge = 1;
    fprintf(1,'Automatic refinement of the corner location after manual mouse click\n');
    disp('Window size for corner finder (wintx and winty):');
    fprintf(1,'wintx ([] = 5) = ');
    wintx = input('');
    if isempty(wintx), wintx = 5; end;
    wintx = round(wintx);
    fprintf(1,'winty ([] = 5) = ');
    winty = input('');
    if isempty(winty), winty = 5; end;
    winty = round(winty);
    
    fprintf(1,'Window size = %dx%d\n',2*wintx+1,2*winty+1);
else
    q_converge = 0;
    fprintf(1,'No attempt to refine the corner location after manual mouse click\n');
end;

if exist(['wintx_' num2str(kk)]),
    eval(['wintxkk = wintx_' num2str(kk) ';']);
    if ~isempty(wintxkk) & ~isnan(wintxkk),
        eval(['wintx = wintx_' num2str(kk) ';']);
        eval(['winty = winty_' num2str(kk) ';']);
    end;
end;

fprintf(1,'Using (wintx,winty)=(%d,%d) - Window size = %dx%d      (Note: To reset the window size, run script clearwin)\n',wintx,winty,2*wintx+1,2*winty+1);

if (exist('map')~=1), map = gray(256); end;


%disp('WARNING!!! Do not forget to change dX_default and dY_default in click_calib.m!!!')

if exist('dX'),
    dX_default = dX;
end;

if exist('dY'),
    dY_default = dY;
end;

if exist('n_sq_x'),
    n_sq_x_default = n_sq_x;
end;

if exist('n_sq_y'),
    n_sq_y_default = n_sq_y;
end;


if ~exist('dX_default')|~exist('dY_default');
    
    % Setup of JY - 3D calibration rig at Intel (new at Intel) - use units in mm to match Zhang
    dX_default = 50;
    dY_default = 50;
     
      
end;

if (exist('dX')~=1)||(exist('dY')~=1)  % This question is now asked only once
    % Enter the size of each square   
    fprintf(1,['Size dX of each square along the X direction ([]=' num2str(dX_default) 'mm) = ']);
    dX = input('');
    fprintf(1,['Size dX of each square along the Y direction ([]=' num2str(dY_default) 'mm) = ']);
    dY = input('');
    if isempty(dX), dX = dX_default; else dX_default = dX; end;
    if isempty(dY), dY = dY_default; else dY_default = dY; end;
else
    fprintf(1,['Size of each square along the X direction: dX=' num2str(dX) 'mm\n']);
    fprintf(1,['Size of each square along the Y direction: dY=' num2str(dY) 'mm   (Note: To reset the size of the squares, clear the variables dX and dY)\n']);
end;


fprintf(1,'Number(s) of image(s) to process ([] = all images) =');
ima_numbers = input('');

if ~exist('n_sq_x_default')|~exist('n_sq_y_default'),
    n_sq_x_default = 10;
    n_sq_y_default = 10;
end;


if isempty(ima_numbers),
    ima_proc = 1:n_ima;
else
    ima_proc = ima_numbers;
end;

for kk = ima_proc,
    
    
    if ~type_numbering,   
        number_ext =  num2str(image_numbers(kk));
    else
        number_ext = sprintf(['%.' num2str(N_slots) 'd'],image_numbers(kk));
    end;
    
    ima_name = [calib_name  number_ext '.' format_image];
    
    
    if exist(ima_name),
        
        fprintf(1,'\nProcessing image %d...\n',kk);

        fprintf(1,'Loading image %s...\n',ima_name);
        
        if format_image(1) == 'p',
            if format_image(2) == 'p',
                I = double(loadppm(ima_name));
            else
                I = double(loadpgm(ima_name));
            end;
        else
            if format_image(1) == 'r',
                I = readras(ima_name);
            else
                I = double(imread(ima_name));
            end;
        end;
        
        
        if size(I,3)>1,
            I = 0.299 * I(:,:,1) + 0.5870 * I(:,:,2) + 0.114 * I(:,:,3);
        end;
        
        eval(['I_' num2str(kk) ' = I;']);
        
        [ny,nx,junk] = size(I);
        Wcal = nx; % to avoid errors later
        Hcal = ny; % to avoid errors later
        
        click_ima_calib_fisheye_no_read;
        %manual_corner_extraction_my;
        active_images(kk) = 1;
        
    else
        eval(['dX_' num2str(kk) ' = NaN;']);
        eval(['dY_' num2str(kk) ' = NaN;']);  
        
        eval(['wintx_' num2str(kk) ' = NaN;']);
        eval(['winty_' num2str(kk) ' = NaN;']);
        
        eval(['x_' num2str(kk) ' = NaN*ones(2,1);']);
        eval(['X_' num2str(kk) ' = NaN*ones(3,1);']);
        
        eval(['n_sq_x_' num2str(kk) ' = NaN;']);
        eval(['n_sq_y_' num2str(kk) ' = NaN;']);
    end;
end;


check_active_images;

% Fix potential non-existing variables:

for kk = 1:n_ima,
    if ~exist(['x_' num2str(kk)]),
        eval(['dX_' num2str(kk) ' = NaN;']);
        eval(['dY_' num2str(kk) ' = NaN;']);  
        
        eval(['x_' num2str(kk) ' = NaN*ones(2,1);']);
        eval(['X_' num2str(kk) ' = NaN*ones(3,1);']);
        
        eval(['n_sq_x_' num2str(kk) ' = NaN;']);
        eval(['n_sq_y_' num2str(kk) ' = NaN;']);
    end;
    
    if ~exist(['wintx_' num2str(kk)]) | ~exist(['winty_' num2str(kk)]),
        
        eval(['wintx_' num2str(kk) ' = NaN;']);
        eval(['winty_' num2str(kk) ' = NaN;']);
        
    end;
end;


string_save = 'save calib_data active_images ind_active wintx winty n_ima type_numbering N_slots first_num image_numbers format_image calib_name Hcal Wcal nx ny map dX_default dY_default dX dY';

for kk = 1:n_ima,
    string_save = [string_save ' X_' num2str(kk) ' x_' num2str(kk) ' n_sq_x_' num2str(kk) ' n_sq_y_' num2str(kk) ' wintx_' num2str(kk) ' winty_' num2str(kk) ' dX_' num2str(kk) ' dY_' num2str(kk)];
end;

eval(string_save);

disp('done');


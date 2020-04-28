
function [y_clean] = ck_rc3_linear_regression_vector(y_noisy,noise)
 %%%%%PERFORM ACTUAL LINEAR REGRESSION REMOVAL OF NOISE COMPONENTS
 %%%%%INPUT: Y=2D Data (1 x col - timepoints)
 %%%%%       noise = noise components (row - components, col - timepoints
 %%%%%OUTPUT: data with removed linear compoents using lin regression 
 %%%%%        mean is not removed

 
 %%%%%%%%%%%%%%%%%%%%
 %%%% Pruefe und korrigiere ggf. die Vektorausrichtung
    dim2    = size(y_noisy,2);
    dim1    = size(y_noisy,1);
    if dim1>dim2
        y_noisy = y_noisy';
        dim2    = size(y_noisy,2);
        dim1    = size(y_noisy,1);
    end
    if dim1~= 1
        fprintf('input vector has to be of size 1 x n\n');
        error('');
    end
    
    ndim1 = size(noise,1);
    ndim2 = size(noise,2);
    
    if ndim1~=dim2
        noise = noise';
        ndim1 = size(noise,1);
        ndim2 = size(noise,2);
    end
    if ndim1~=dim2
        fprintf('noise vector has is not of correct size\n');
        error('');
    end
    %%%% Pruefung abgeschlossen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    assignin('base','y_noisyx',y_noisy)
    assignin('base','noisex',noise)
    
    
    %%% build the X and add the b0
    X       = [];
    X = [ones(dim2,1)];
    X = [X noise];
     assignin('base','Xx',X)
   

    %%%% run through all voxels and perform  linear regression

    beta = X\y_noisy';
%    beta(1,:) = 1; %%% move the mw
    beta(1,:) = 0; %%%do not move the mw
    
            yrp = X * beta;
            y_clean=y_noisy - yrp';
    

end
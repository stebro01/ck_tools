function [x y z stat] = ck_rc3_mni2nifti_real_with_output(V, xi, yi, zi,stat)
% wenn stat == 1 dann wird ein output generiert
xyz = [xi yi zi 1] * inv(V(1,1).mat') ;
x = round(xyz(1));
y = round(xyz(2));
z = round(xyz(3));
if abs(x-xyz(1))>0.01 && stat
    fprintf('fehler in ck_rc2_mni2nifti_real\n');
    fprintf(' x = %.1f   y= %.1f  z = %.1f \n',x,y,z);
    fprintf(' xyz(1) = %.1f   xyz(2)= %.1f  xyz(3)in = %.1f \n',xyz(1),xyz(2),xyz(3));
    fprintf(' xin = %.1f   yin= %.1f  zin = %.1f \n',xi,yi,zi);
    assignin('base','Vi',V);
    stat =0;
end
if abs(y-xyz(2))>0.01 && stat
    fprintf('fehler in ck_rc2_mni2nifti_real\n');
    fprintf(' x = %.1f   y= %.1f  z = %.1f \n',x,y,z);
    fprintf(' xin = %.1f   yin= %.1f  zin = %.1f \n',xi,yi,zi);
    assignin('base','Vi',V);
    stat =0;
end
if abs(z-xyz(3))>0.01 && stat
    fprintf('fehler in ck_rc2_mni2nifti_real\n');
    fprintf(' x = %.1f   y= %.1f  z = %.1f \n',x,y,z);
    fprintf(' xin = %.1f   yin= %.1f  zin = %.1f \n',xi,yi,zi);
    stat = 0;
end

% Text von John Ashburner
% mm  = [20 30 -10 ; 40 -20 30]; % example co-ordinates.
% V   = spm_vol(spm_get(1,'*.img'));
% M   = inv(V.mat);
% vox = (M(1:3,1:3)*mm' + repmat(M(1:3,4),1,size(mm,1)))'
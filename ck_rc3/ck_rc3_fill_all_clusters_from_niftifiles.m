function stat = ck_rc3_fill_all_clusters_from_niftifiles(clusterdir,outdir,Vm,Pm,options)
% Füllt die im Verzeichnis stehenden Clusterstrukturen mit den
% Daten aus den Files

%options
%Pm
%Vm
if options.extract_con1
    c11val=options.extract_con1_g1_val;
    c12val=options.extract_con1_g2_val;
    c13val=options.extract_con1_g3_val;
end

if options.extract_con2
    c21val=options.extract_con2_g1_val;
    c22val=options.extract_con2_g2_val;
    c23val=options.extract_con2_g3_val;
end

do_append = 0; %ADDED by SB280518
if isfield(options, 'method_str')
   switch lower(options.method_str)
       case 'append'
           do_append = 1;
   end
end

stat = 0;
load(fullfile(clusterdir,'parameter.mat'));
if do_append == 0
   options.cluster_fn = fullfile(clusterdir,'Cluster.mat');
end
load(options.cluster_fn); fprintf('|Loading Cluster: %s\n', options.cluster_fn);

h =  waitbar(0,'Please wait... filling Clusters');
% Schleife ueber die Probanden
for i=1:length(Vm)
    waitbar(i/length(Vm));
    cluster_save_name=fullfile(outdir,['Cluster_' Pm{i} '.mat']);
    %data_save_name = fullfile(outdir,['Data_' Pm{i,1} '.mat']);
    
    V=spm_vol(Vm{i});
    assignin('base','Vm',Vm);
    assignin('base','V',V);
    [Y,Data.XYZ]=spm_read_vols(V);
    assignin('base','Y',Y);
    assignin('base','XYZ',Data.XYZ);
    % falls option gesetzt dann grenze den Cluster entsprechend den
    % Bedingungen ein
    if options.extract_con1 || options.extract_con2
        [p,f]=fileparts(Vm{i}(1,:));
        fc1 = fullfile(p,options.extract_con_g1_name);
        fc2 = fullfile(p,options.extract_con_g2_name);
        fc3 = fullfile(p,options.extract_con_g3_name);
        Vc1 = spm_vol(fc1);
        Vc2 = spm_vol(fc2);
        Vc3 = spm_vol(fc3);
        [Ycc1,~]=spm_read_vols(Vc1);
        [Ycc2,~]=spm_read_vols(Vc2);
        [Ycc3,~]=spm_read_vols(Vc3);
    end
    % Schleife ueber die Cluster
    for c=1:size(Cluster,1) %#ok<NODEF,USENS>
        jc = 0;
        % zaehler fuer den conditional Cluster
        jcc1 = 0;
        jcc2 = 0;
        num_vox = size(Cluster{c,1}.XYZ,2);
        Yall = zeros(size(Y,4),num_vox); %#ok<AGROW>
        Yc1 = zeros(size(Y,4),num_vox); %#ok<AGROW>
        Yc2 = zeros(size(Y,4),num_vox); %#ok<AGROW>
        stat = 1;
        % Schleife ueber die Voxel
        for j=1:num_vox
            
            x_mni = Cluster{c,1}.XYZ(1,j);
            y_mni = Cluster{c,1}.XYZ(2,j);
            z_mni = Cluster{c,1}.XYZ(3,j);
            %               [x,y,z,stat] = ck_rc2_mni2nifti_real(V(1,1),x_mni,y_mni,z_mni);
            [x,y,z,stat] = ck_rc3_mni2nifti_real_with_output(V(1,1),x_mni,y_mni,z_mni,stat);
            
            if stat ==0
                stat = 1;
                fprintf('Fehler bei Cluster nummer %.0f von Proband %s\n',c,Vm{i,1}(1,:));
            end
            
            % kopieren der Zeitreihe
            if x > 0 && y > 0 && z > 0
                jc = jc +1;
                Yall(:,jc) = squeeze(Y(x,y,z,:)); %#ok<AGROW>
                if options.extract_con1
                    if Ycc1(x,y,z)>c11val && Ycc2(x,y,z)<c12val && Ycc3(x,y,z)<c13val
                        jcc1 = jcc1 + 1;
                        Yc1(:,jcc1) =  Yall(:,jc);
                    end
                end
                if options.extract_con2
                    if Ycc1(x,y,z)>c21val && Ycc2(x,y,z)<c22val && Ycc3(x,y,z)<c23val
                        jcc2 = jcc2 + 1;
                        Yc2(:,jcc2) =  Yall(:,jc);
                    end
                end
            else
                fprintf('neg indices eigentlich nicht moeglich x = %.1f y=%.1f z=%.1f\n');
            end
        end
        % alle Voxel des clusters
        Cluster{c,1}.Y = Yall(:,1:jc);
        Cluster{c,1}.Y_mean = mean(Cluster{c,1}.Y, 2, 'omitnan'); %#ok<AGROW> %OMITNAN: added by SB 20.06.18
        if sum(Cluster{c,1}.Y_mean) == 0, disp('please check the Y_mean, sb_092342390!'), return, end
        if options.extract_con1
            % nur Voxel welche die in options angegebenen Bedinungen
            % erfuellen
            Cluster{c,1}.Yc1 = Yc1(:,1:jcc1);
            Cluster{c,1}.Yc1_mean = mean(Cluster{c,1}.Yc1, 2, 'omitnan'); %#ok<AGROW> %OMITNAN: added by SB 20.06.18
            %fprintf('conditional .... %.0f / %.0f / %.0f \n',num_vox,jc,jcc)
        end
        if options.extract_con2
            % nur Voxel welche die in options angegebenen Bedinungen
            % erfuellen
            Cluster{c,1}.Yc2 = Yc2(:,1:jcc2);
            Cluster{c,1}.Yc2_mean = mean(Cluster{c,1}.Yc2, 2, 'omitnan'); %#ok<AGROW>%OMITNAN: added by SB 20.06.18
            %fprintf('conditional .... %.0f / %.0f / %.0f \n',num_vox,jc,jcc)
        end
    end
    fprintf('speichere nun unter folgendem Filename: %s\n',cluster_save_name);
    
    %added by sb280518: try to append Cluster
    if do_append
        if exist(cluster_save_name)
            fprintf('\t|%s found, appending new data\n', cluster_save_name);
           Cluster_old = load(cluster_save_name);
        else
            fprintf('\t|%s not found, creating Cluster new\n', cluster_save_name);
            Cluster_old.Cluster = [];
        end
        
        Cluster_old.Cluster = [Cluster_old.Cluster;Cluster];
        Cluster = Cluster_old.Cluster;
    end %do_append
    %end added sb280518
    
    save(cluster_save_name,'Cluster','-v7.3');
    %save(data_save_name,'Data');
    
    
end
close(h)

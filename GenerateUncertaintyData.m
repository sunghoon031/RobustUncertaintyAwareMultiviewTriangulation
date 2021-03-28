clear all; close all; clc;

load_precomputed_results = true;

if (load_precomputed_results)
    load uncertainty.mat;
    load results_inliers.mat;
else

    tic
    % Camera follows TUM rgb-d dataset parameters
    focal_length = 525;
    img_width = 640;
    img_height = 480;
    K = [focal_length 0 0; 0 focal_length 0; 0 0 1];

    noise_levels = [1:10, 12:2:30];
    dists = [1:9, 10:10:30];

    n_cameras = [2:10, 12:2:30, 35:5:50];
    
    n_simulations = [];
    for i = 1:length(n_cameras)
        
        if (n_cameras(i) <= 5)
            n_simulations = [n_simulations, 10];
            continue;
        end
        
        if (n_cameras(i) <= 20)
            n_simulations = [n_simulations, 5];
            continue;
        end
        
        n_simulations = [n_simulations, 3];
    end
    
    n_simulations = 1000*n_simulations;
    
    n_total_run = length(noise_levels)*length(dists)*sum(n_simulations);

    all_error_3D = nan(1,n_total_run);
    all_error_2D = nan(1,n_total_run);
    all_min_level_of_degeneracy = nan(1,n_total_run);
    all_n_camera = nan(1,n_total_run);

    idx = 0;

    for noise_level = noise_levels

        for dist = dists
            point_w = [0;0;dist];

            for n_camera = fliplr(n_cameras)

                disp(['noise_level = ', num2str(noise_level), ', dist = ', num2str(dist), ', n_camera = ', num2str(n_camera)])


                n_simulation = n_simulations(n_cameras == n_camera);
                for simulation = 1:n_simulation

                    idx = idx + 1;
                    cam_pos = cell(1,n_camera);
                    cam_pos_mat = nan(3,n_camera);
                    R = cell(1,n_camera);
                    t = cell(1,n_camera);
                    uv = cell(1,n_camera);
                    uv_measured = cell(1,n_camera);
                    ray = cell(1,n_camera);
                    AA = cell(1,n_camera);
                    bb = cell(1,n_camera);

                    for i = 1:n_camera

                        % Camera is located randomly inside a unit sphere.
                        cam_pos{i} = rand(3,1)-0.5;
                        cam_pos{i} = rand(1)*0.5*cam_pos{i}/norm(cam_pos{i});

                        % The first two cameras are antipodal.
                        if (i==1)
                            cam_pos{1} = rand(3,1)-0.5;
                            cam_pos{1} = 0.5*cam_pos{1}/norm(cam_pos{1});
                        elseif (i==2)
                            cam_pos{2} = -cam_pos{1};
                        end

                        while(true)
                            R_temp = random_rotation(180);
                            t_temp = -R_temp*cam_pos{i};

                            point_c = R_temp*point_w+t_temp;
                            uv_temp = K*[point_c(1)/point_c(3); point_c(2)/point_c(3);1]; 
                            uv_temp = uv_temp(1:2);
                            if (abs(uv_temp(1)) > img_width/2 || abs(uv_temp(2)) > img_height/2 || point_c(3) < 0)
                                continue;
                            end

                            while(true)
                                noise_dir = rand(2,1)-0.5;
                                noise_dir = noise_dir/norm(noise_dir);
                                uv_measured_temp = uv_temp+normrnd(0, noise_level)*noise_dir;
                                if (abs(uv_measured_temp(1)) > img_width/2 || abs(uv_measured_temp(2)) > img_height/2)
                                    continue;
                                end
                                break;
                            end


                            break;
                        end

                        cam_pos_mat(:,i) = cam_pos{i};
                        R{i} = R_temp;
                        t{i} = t_temp;
                        uv{i} = uv_temp;
                        uv_measured{i} = uv_measured_temp;
                        ray{i} = K\[uv_measured{i};1];
                        ray{i} = ray{i}/norm(ray{i});

                        % Compute AA and bb necesary for computing the Jacobian:
                        r11 = R{i}(1,1); r12 = R{i}(1,2); r13 = R{i}(1,3);
                        r21 = R{i}(2,1); r22 = R{i}(2,2); r23 = R{i}(2,3);
                        r31 = R{i}(3,1); r32 = R{i}(3,2); r33 = R{i}(3,3);
                        tx = t{i}(1); ty = t{i}(2); tz = t{i}(3);
                        fx = K(1,1); skew = K(1,2); fy = K(2,2); 

                        a11 = [...
                            0, ...
                            fx*(r11*r32-r31*r12)+skew*(r21*r32-r31*r22), ...
                            fx*(r11*r33-r31*r13)+skew*(r21*r33-r31*r23), ...
                            fx*(r11*tz-r31*tx)+skew*(r21*tz-r31*ty)];
                        a21 = [...
                            0, ...
                            fy*(r21*r32-r31*r22), ...
                            fy*(r21*r33-r31*r23), ...
                            +fy*(r21*tz-r31*ty)];
                        a12 = [...
                            fx*(r12*r31-r32*r11)+skew*(r22*r31-r32*r21),...
                            0,...
                            fx*(r12*r33-r32*r13)+skew*(r22*r33-r32*r23),...
                            fx*(r12*tz-r32*tx)+skew*(r22*tz-r32*ty)];
                        a22 = [...
                            fy*(r22*r31-r32*r21),...
                            0,...
                            fy*(r22*r33-r32*r23),...
                            fy*(r22*tz-r32*ty)];
                        a13 = [...
                            fx*(r13*r31-r33*r11)+skew*(r23*r31-r33*r21),...
                            fx*(r13*r32-r33*r12)+skew*(r23*r32-r33*r22),...
                            0,...
                            fx*(r13*tz-r33*tx)+skew*(r23*tz-r33*ty)];
                        a23 = [...
                            fy*(r23*r31-r33*r21),...
                            fy*(r23*r32-r33*r22),...
                            0,...
                            fy*(r23*tz-r33*ty)];

                        AA{i} = [a11;a21;a12;a22;a13;a23];
                        bb{i} = [r31, r32, r33, tz];
                    end


        %             figure;
        %             axis equal
        %             hold on
        %             for i=1:n_camera
        %                 if (i<= n_inlier)
        %                     scatter3(cam_pos{i}(1),cam_pos{i}(2),cam_pos{i}(3), 'ko')
        %                     text(cam_pos{i}(1), cam_pos{i}(2), cam_pos{i}(3),num2str(i),'HorizontalAlignment','left','FontSize',20);   
        %                     
        %                     endpoint = 15*R{i}'*ray{i}+cam_pos{i};
        %                     xx = [cam_pos{i}(1), endpoint(1)];
        %                     yy = [cam_pos{i}(2), endpoint(2)];
        %                     zz = [cam_pos{i}(3), endpoint(3)];
        %                     plot3(xx, yy, zz, 'Color', [0 0 0 0.5]);
        %                 else
        %                     scatter3(cam_pos{i}(1),cam_pos{i}(2),cam_pos{i}(3), 'ro')
        %                     
        %                     endpoint = R{i}'*ray{i}+cam_pos{i};
        %                     xx = [cam_pos{i}(1), endpoint(1)];
        %                     yy = [cam_pos{i}(2), endpoint(2)];
        %                     zz = [cam_pos{i}(3), endpoint(3)];
        %                     plot3(xx, yy, zz, 'Color', [1 0 0 0.5]);
        %                 end
        %             end
        %             scatter3(point_w(1), point_w(2), point_w(3), 200, 'kx')
        %             view(360,0)


                    %% DLT + GN for all inliers

                    % Prepare the matrices.
                    M0 = nan(1,n_camera);
                    M1 = nan(1,n_camera);    
                    M2 = nan(4,n_camera);
                    M3 = nan(4,n_camera);
                    M4 = nan(4,n_camera);
                    ray_w = cell(1,n_camera);
                    for i = 1:n_camera
                        ray_w{i} = K\[uv_measured{i};1];
                        ray_w{i} = ray_w{i}/norm(ray_w{i});
                        ray_w{i} = R{i}'*ray_w{i};

                        M0(1,i) = uv_measured{i}(1)-K(1,3);
                        M1(1,i) = uv_measured{i}(2)-K(2,3);

                        P = [R{i}, t{i}];
                        M2(:,i) = K(1,1)*P(1,:)'+K(1,2)*P(2,:)';
                        M3(:,i) = K(2,1)*P(1,:)'+K(2,2)*P(2,:)';
                        M4(:,i) = P(3,:)';
                    end

                    % Do DLT.
                    A = zeros(2*n_camera, 4);
                    for i = 1:n_camera
                        Q = K*[R{i}, t{i}];
                        u = uv_measured{i}(1);
                        v = uv_measured{i}(2);
                        A(2*i-1:2*i, :) = [u*Q(3,:)-Q(1,:); v*Q(3,:)-Q(2,:)];
                    end
                    [~,S_,V_] = svd(A);
                    [~, min_idx] = min(diag(S_));
                    x_dlt = V_(:, min_idx);
                    x_dlt = x_dlt/x_dlt(end);
                    x_dlt = x_dlt(1:3);

                    % Do Gauss-Newton Optimization.
                    mean_reproj_error = 0;
                    x_gn = x_dlt;
                    for it = 1:10
                        N4 = [x_gn;1]'*M4;
                        M5 = ([x_gn;1]'*M2)./N4 - M0;
                        M6 = ([x_gn;1]'*M3)./N4 - M1;
                        M5 = M5(1:n_camera);
                        M6 = M6(1:n_camera);
                        residuals = [M5;M6];
                        squared_reproj_error_all = M5.*M5+M6.*M6;

                        mean_reproj_error_prev = mean_reproj_error;
                        mean_reproj_error = mean(sqrt(squared_reproj_error_all));

                        if (it > 2 && abs(mean_reproj_error-mean_reproj_error_prev) < 0.1)
                            break;
                        end

                        J_est = zeros(2*n_camera,3);
                        for i = 1:n_camera
                            J_temp = AA{i}*[x_gn;1]/(bb{i}*[x_gn;1])^2;
                            J_est(2*i-1:2*i, :) = reshape(J_temp, [2,3]); 
                        end
                        update = -pinv(J_est)*residuals(:);
                        x_gn = x_gn+update;
                    end

                    % Evaluate the performance metrics.
                    N4 = [x_gn;1]'*M4;
                    M5 = ([x_gn;1]'*M2)./N4 - M0;
                    M6 = ([x_gn;1]'*M3)./N4 - M1;
                    squared_reproj_error_all = M5.*M5+M6.*M6;
                    squared_reproj_error_inliers = squared_reproj_error_all(1:n_camera);


                    all_error_3D(idx) = norm(point_w-x_gn);
                    all_error_2D(idx) = mean(sqrt(squared_reproj_error_inliers));
                    all_min_level_of_degeneracy(idx) = GetMinimumLevelOfDegeneracy(cam_pos, ray_w, 1:n_camera);
                    all_n_camera(idx) = n_camera;

                end

            end
        end
    end
    
    angle_bin_edges = fliplr([0:19; 1:20]);
    noise_bin_edges = [0:19; 1:20];
    cell_colormap_results = cell(size(angle_bin_edges,2), size(noise_bin_edges,2), length(n_cameras));

    for i = 1:length(all_error_3D)
        angle_of_degeneracy = acosd(all_min_level_of_degeneracy(i));
        n_camera = all_n_camera(i);
        error_3D = all_error_3D(i);
        error_2D = all_error_2D(i);
        for j = 1:size(angle_bin_edges,2)
            if (angle_of_degeneracy < angle_bin_edges(1,j) || angle_of_degeneracy >= angle_bin_edges(2,j))
                continue;
            end
            for k = 1:size(noise_bin_edges,2)
                if (error_2D < noise_bin_edges(1,k) || error_2D >= noise_bin_edges(2,k))
                    continue;
                end
                for l = 1:length(n_cameras)
                    if (n_cameras(l) ~= n_camera)
                        continue;
                    end
                    cell_colormap_results{j,k,l} = [cell_colormap_results{j,k,l}, error_3D];
                end
            end
        end
    end

    mat_colormap_results = zeros(size(cell_colormap_results));
    for j = 1:size(mat_colormap_results,1)
        for k = 1:size(mat_colormap_results,2)
            for l = 1:size(mat_colormap_results,3)
                mat_colormap_results(j,k,l) = min(1,sqrt(mean(cell_colormap_results{j,k,l}.^2)));
            end
        end
    end
    
    mat_colormap_results_smoothed = MonotoneSmooth(mat_colormap_results, 0.01, 1, 1, 0);

    uncertainty.angle_i = angle_bin_edges;
    uncertainty.reproj_j = noise_bin_edges;
    uncertainty.camera_k = n_cameras;
    uncertainty.result = mat_colormap_results_smoothed;

    save uncertainty.mat uncertainty
    save results_inliers.mat
    time_elsaped = toc;
    disp(['Took ', num2str(time_elsaped/3600), ' hours.'])
    
end



%out = TrilinearInterpolation(8.5, 3.5, 10, uncertainty)





%% Plot 

figure;
for i = 1:length(n_cameras)
    subplot(6,8,2*i-1)
    imagesc(mat_colormap_results(:,:,i))
    colormap(jet)
    caxis([0 1])
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    title([num2str(n_cameras(i)), ' cams'])
    
    subplot(6,8,2*i)
    imagesc(mat_colormap_results_smoothed(:,:,i))
    colormap(jet)
    caxis([0 1])
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    title([num2str(n_cameras(i)), ' cams (smoothed)'])
end


figure;
plot_idx = 0;
for i = 1:length(n_cameras)
    n_camera = n_cameras(i);
    if (ismember(n_cameras(i), [2, 3, 4, 5, 10, 20,30,50]))
        plot_idx = plot_idx + 1;
        subplot(2,4,plot_idx)
        imagesc(mat_colormap_results_smoothed(:,:,i))
        axis equal
        colormap(jet)
        colorbar
        caxis([0 1])
        xticks(0.5:5:20.5)
        xticklabels({'0', '5', '10', '15', '20'})
        xlabel('Mean 2D error (pix)')
        yticks(0.5:5:20.5)
        yticklabels({'20','15', '10', '5','0'})
        ylabel('Angle of Degeneracy (deg)')
        title([num2str(n_cameras(i)), ' cams'])
    end
end

% for i = 1:length(n_cameras)
%     disp(['n_camera = ', num2str(n_cameras(i))])
%     latex_table_raw = latex(vpa(sym(mat_colormap_results(:,:,i)),2))
%     latex_table_smoothed = latex(vpa(sym(mat_colormap_results_smoothed(:,:,i)),2))
% end


    %%
function out = random_rotation(max_angle_deg)

    axis = rand(3,1)-0.5;
    axis = axis/norm(axis);
    angle = rand*max_angle_deg/180*pi;
    rotvec = angle*axis;
    out = rotation_from_rotvec(rotvec);

end

function R = rotation_from_rotvec(in)
    angle = norm(in);
    
    if (angle==0)
        R = eye(3);
    else
        unit_axis = in/angle;
        so3 = SkewSymmetricMatrix(unit_axis);
        R = eye(3)+so3*sin(angle)+so3^2*(1-cos(angle));
    end
end


function out = SkewSymmetricMatrix(in)
    out=[0 -in(3) in(2) ; in(3) 0 -in(1) ; -in(2) in(1) 0 ];
end


function mLoD = GetMinimumLevelOfDegeneracy(cam_pos, ray_w, idx)
    
    idx = reshape(idx, [1, numel(idx)]);
    
    mLoD = inf;
    for i = idx
        for j = idx
            if (j<=i)
                continue;
            end
            t_ij = cam_pos{i}-cam_pos{j};
            t_ij = t_ij/norm(t_ij);
            f_i = ray_w{i};
            f_j = ray_w{j};
            
            LoD = abs(f_i'*f_j);
%             LoD = max(abs([t_ij'*f_i, t_ij'*f_j, f_i'*f_j]));
            if (LoD < mLoD)
                mLoD = LoD;
            end
        end
    end
    
end


function out = MonotoneSmooth(in, thr, x_increasing, y_increasing, z_increasing)

    out = in;
    for it = 1:100000
        out_prev = out;

        for i = 1:size(out_prev,1)
            for j = 1:size(out_prev,2)
                for k = 1:size(out_prev,3)
                    candidates = [];

                    if (x_increasing)
                        if (i-1 >= 1)
                            if (out_prev(i-1,j,k) <= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i-1,j,k))];
                            end
                        end
                        if (i+1 <= size(out_prev,1))
                            if (out_prev(i+1,j,k) >= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i+1,j,k))];
                            end
                        end
                    else
                        if (i-1 >= 1)
                            if (out_prev(i-1,j,k) >= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i-1,j,k))];
                            end
                        end
                        if (i+1 <= size(out_prev,1))
                            if (out_prev(i+1,j,k) <= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i+1,j,k))];
                            end
                        end
                    end

                    if (y_increasing)
                        if (j-1 >= 1)
                            if (out_prev(i,j-1,k) <= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j-1,k))];
                            end
                        end
                        if (j+1 <= size(out_prev,2))
                            if (out_prev(i,j+1,k) >= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j+1,k))];
                            end
                        end
                    else
                        if (j-1 >= 1)
                            if (out_prev(i,j-1,k) >= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j-1,k))];
                            end
                        end
                        if (j+1 <= size(out_prev,2))
                            if (out_prev(i,j+1,k) <= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j+1,k))];
                            end
                        end
                    end
                    
                    if (z_increasing)
                        if (k-1 >= 1)
                            if (out_prev(i,j,k-1) <= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j,k-1))];
                            end
                        end
                        if (k+1 <= size(out_prev,3))
                            if (out_prev(i,j,k+1) >= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j,k+1))];
                            end
                        end
                    else
                        if (k-1 >= 1)
                            if (out_prev(i,j,k-1) >= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j,k-1))];
                            end
                        end
                        if (k+1 <= size(out_prev,3))
                            if (out_prev(i,j,k+1) <= out_prev(i,j,k))
                                candidates = [candidates, out_prev(i,j,k)];
                            else
                                candidates = [candidates, 0.5*(out_prev(i,j,k)+out_prev(i,j,k+1))];
                            end
                        end
                    end

                    out(i,j,k) = mean(candidates);
                end
            end
        end

             
        % Check if monotonicity is enforced for every value
        
        completed = true;
        for i = 1:size(out,1)
            if (~completed)
                break;
            end
            for j = 1:size(out,2)
                if (~completed)
                    break;
                end
                for k = 1:size(out,3)
                    if (~completed)
                        break;
                    end
                    
                    if (x_increasing)
                        if (out(max(1,i-1),j,k) > out(i,j,k)+thr)
                            completed = false;
                            break;
                        end
                        if (out(min(size(out,1),i+1),j,k) < out(i,j,k)-thr)
                            completed = false;
                            break;
                        end
                    else
                        if (out(max(1,i-1),j,k) < out(i,j,k)-thr)
                            completed = false;
                            break;
                        end
                        if (out(min(size(out,1),i+1),j,k) > out(i,j,k)+thr)
                            completed = false;
                            break;
                        end
                    end
                    
                    if (y_increasing)
                        if (out(i,max(1,j-1),k) > out(i,j,k)+thr)
                            completed = false;
                            break;
                        end
                        if (out(i,min(size(out,2),j+1),k) < out(i,j,k)-thr)
                            completed = false;
                            break;
                        end
                    else
                        if (out(i,max(1,j-1),k) < out(i,j,k)-thr)
                            completed = false;
                            break;
                        end
                        if (out(i,min(size(out,2),j+1),k) > out(i,j,k)+thr)
                            completed = false;
                            break;
                        end
                    end
                    
                    if (z_increasing)
                        if (out(i,j,max(1,k-1)) > out(i,j,k)+thr)
                            completed = false;
                            break;
                        end
                        if (out(i,j,min(size(out,3),k+1)) < out(i,j,k)-thr)
                            completed = false;
                            break;
                        end
                    else
                        if (out(i,j,max(1,k-1)) < out(i,j,k)-thr)
                            completed = false;
                            break;
                        end
                        if (out(i,j,min(size(out,3),k+1)) > out(i,j,k)+thr)
                            completed = false;
                            break;
                        end
                    end
                    
                end
            end
        end
        
        
        
        
        if (completed)
            it
            break;
        end
        
        
        

    end

end


function [out] = TrilinearInterpolation(angle, reproj, n_cam, uncertainty)

    if (reproj > max(uncertainty.reproj_j(2,:)))
        out = 1;
        return;
    end
    
    angle_i = mean(uncertainty.angle_i);
    reproj_j = mean(uncertainty.reproj_j);
    camera_k = uncertainty.camera_k;

    [angle_lb, i_lb] = min(angle_i);
    [angle_ub, i_ub] = max(angle_i);
    for i = 1:length(angle_i)
        angle_candidate = angle_i(i);
        if (angle-angle_candidate >= 0 && angle_candidate > angle_lb)
            angle_lb = angle_candidate;
            i_lb = i;
        end
        if (angle-angle_candidate <= 0 && angle_candidate < angle_ub)
            angle_ub = angle_candidate;
            i_ub = i;
        end
    end

    [reproj_lb, j_lb] = min(reproj_j);
    [reproj_ub, j_ub] = max(reproj_j);
    for j = 1:length(reproj_j)
        reproj_candidate = reproj_j(j);
        if (reproj-reproj_candidate >= 0 && reproj_candidate > reproj_lb)
            reproj_lb = reproj_candidate;
            j_lb = j;
        end
        if (reproj-reproj_candidate <= 0 && reproj_candidate < reproj_ub)
            reproj_ub = reproj_candidate;
            j_ub = j;
        end
    end

    [n_cam_lb, k_lb] = min(camera_k);
    [n_cam_ub, k_ub] = max(camera_k);
    for k = 1:length(camera_k)
        n_cam_candidate = camera_k(k);
        if (n_cam-n_cam_candidate >= 0 && n_cam_candidate > n_cam_lb)
            n_cam_lb = n_cam_candidate;
            k_lb = k;
        end
        if (n_cam-n_cam_candidate <= 0 && n_cam_candidate < n_cam_ub)
            n_cam_ub = n_cam_candidate;
            k_ub = k;
        end
    end


    if (angle_ub == angle_lb)
        x_d = 0;
    else
        x_d = (angle-angle_lb)/(angle_ub-angle_lb);
    end
    if (reproj_ub == reproj_lb)
        y_d = 0;
    else
        y_d = (reproj-reproj_lb)/(reproj_ub-reproj_lb);
    end
    if (n_cam_ub == n_cam_lb)
        z_d = 0;
    else
        z_d = (n_cam-n_cam_lb)/(n_cam_ub-n_cam_lb);
    end

    c_000 = uncertainty.result(i_lb,j_lb,k_lb);
    c_001 = uncertainty.result(i_lb,j_lb,k_ub);
    c_010 = uncertainty.result(i_lb,j_ub,k_lb);
    c_011 = uncertainty.result(i_lb,j_ub,k_ub);
    c_100 = uncertainty.result(i_ub,j_lb,k_lb);
    c_101 = uncertainty.result(i_ub,j_lb,k_ub);
    c_110 = uncertainty.result(i_ub,j_ub,k_lb);
    c_111 = uncertainty.result(i_ub,j_ub,k_ub);

    c_00 = c_000*(1-x_d)+c_100*x_d;
    c_01 = c_001*(1-x_d)+c_101*x_d;
    c_10 = c_010*(1-x_d)+c_110*x_d;
    c_11 = c_011*(1-x_d)+c_111*x_d;

    c_0 = c_00*(1-y_d)+c_10*y_d;
    c_1 = c_01*(1-y_d)+c_11*y_d;

    out = c_0*(1-z_d)+c_1*z_d;
    
%     disp(['angle = ', num2str(angle), ', lb = ', num2str(angle_lb), ', ub = ', num2str(angle_ub), ' angle(i_lb) = ', num2str(angle_i(i_lb)), ', angle(i_ub) = ', num2str(angle_i(i_ub))])
%     disp(['reproj = ', num2str(reproj), ', lb = ', num2str(reproj_lb), ', ub = ', num2str(reproj_ub), ' reproj(j_lb) = ', num2str(reproj_j(j_lb)), ', reproj(j_ub) = ', num2str(reproj_j(j_ub))])
%     disp(['n_cam = ', num2str(n_cam), ', lb = ', num2str(n_cam_lb), ', ub = ', num2str(n_cam_ub), ' camera_k(k_lb) = ', num2str(camera_k(k_lb)), ', camera_k(k_ub) = ', num2str(camera_k(k_ub))])


end
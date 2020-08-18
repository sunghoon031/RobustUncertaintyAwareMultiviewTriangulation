clc;

%% Load the 3D uncertainty grid
load uncertainty.mat;

%% World point to be triangulated
dist = 2+rand(1)*8; % between 2 and 10 unit
point_w = [0;0;dist];

%% Set up the cameras (50 cameras with inlying measurements + 50 cameras with outlying measurements)
focal_length = 525;
img_width = 640;
img_height = 480;

noise_level = 3; % pixel (Gaussian noise level)
outlier_min_noise = 10; % pixel (Outliers contain noise greater than this) 

n_inliers = 50;
n_outliers = 50;
n_cameras = n_inliers + n_outliers;

[R, t, K, uv_measured] = SimulateCameras(point_w, n_inliers, n_cameras, img_width, img_height, focal_length, noise_level, outlier_min_noise);

%% Triangulate the point using our method.

% Precomputation indepndent of the point
tic;
[cam_pos, P, A, M3, M4, M5] = Precompute(n_cameras, K, R, t);
time_precompute = toc;

% Gauss-Newton optimization:
tic;
[x_gn, uncertainty_3D_gn] = RunMultiviewTriangulation('GN', uncertainty, cam_pos, R,P,K,uv_measured,A, M3, M4, M5);
time_gn = toc;

% DLT optimization:
tic;
[x_dlt, uncertainty_3D_dlt] = RunMultiviewTriangulation('DLT', uncertainty, cam_pos, R,P,K,uv_measured,A, M3, M4, M5);
time_dlt = toc;


%% Evaluate the performance

error_2D_gn = GetMean2DErrorOfTrueInliers(x_gn, n_cameras, n_inliers, K, P, uv_measured);
error_2D_dlt = GetMean2DErrorOfTrueInliers(x_dlt, n_cameras, n_inliers, K, P, uv_measured);

error_3D_gn = norm(x_gn-point_w);
error_3D_dlt = norm(x_dlt-point_w);

disp(['Point distance = ', num2str(dist), ' unit (1 unit = geometric span of the cameras).'])
disp(['Number of cameras = ', num2str(n_cameras), ', Outlier ratio = ', num2str(n_outliers/n_cameras*100), '%.'])
disp(['Precomputation took ', num2str(time_precompute*1000), ' ms.'])
disp(['GN optimization:  2D error = ', num2str(error_2D_gn), ' pix, 3D error = ', num2str(error_3D_gn), ' unit, uncertainty = ', num2str(uncertainty_3D_gn), ', took ', num2str(time_gn*1000), ' ms.'])
disp(['DLT optimization: 2D error = ', num2str(error_2D_dlt), ' pix, 3D error = ', num2str(error_3D_dlt), ' unit, uncertainty = ', num2str(uncertainty_3D_dlt), ', took ', num2str(time_dlt*1000), ' ms.'])







%% Function definitions

function R = RandomRotation(max_angle_rad)

    unit_axis = rand(3,1)-0.5;
    unit_axis = unit_axis/norm(unit_axis);
    angle = rand*max_angle_rad;
    R = RotationFromUnitAxisAngle(unit_axis, angle);

end

function R = RotationFromUnitAxisAngle(unit_axis, angle)
    
    if (angle==0)
        R = eye(3);
    else
        so3 = SkewSymmetricMatrix(unit_axis);
        R = eye(3)+so3*sin(angle)+so3^2*(1-cos(angle));
    end
end


function out = SkewSymmetricMatrix(in)
    out=[0 -in(3) in(2) ; in(3) 0 -in(1) ; -in(2) in(1) 0 ];
end

function [R, t, K, uv_measured] = SimulateCameras(point_w, n_inlier, n_cameras, img_width, img_height, focal_length, noise_level, outlier_min_noise)

    K_cam = [focal_length 0 0; 0 focal_length 0; 0 0 1];
    
    
    R = cell(1,n_cameras);
    t = cell(1,n_cameras);
    K = cell(1,n_cameras);
    uv_measured = cell(1,n_cameras);

    for i = 1:n_cameras

        K{i} = K_cam;

        % Camera is located randomly inside a unit sphere.
        cam_pos = rand(3,1)-0.5;
        cam_pos = rand(1)*0.5*cam_pos/norm(cam_pos);

        % The first two cameras are antipodal.
        if (i==1)
            cam_pos = rand(3,1)-0.5;
            cam_pos = 0.5*cam_pos/norm(cam_pos);
        elseif (i==2)
            cam_pos = -cam_pos;
        end

        while(true)
            R_temp = RandomRotation(pi);
            t_temp = -R_temp*cam_pos;

            point_c = R_temp*point_w+t_temp;
            uv = K{i}*[point_c(1)/point_c(3); point_c(2)/point_c(3);1]; 
            uv = uv(1:2);
            if (abs(uv(1)) > img_width/2 || abs(uv(2)) > img_height/2 || point_c(3) < 0)
                continue;
            end

            if (i <= n_inlier)
                while(true)
                    noise_dir = rand(2,1)-0.5;
                    noise_dir = noise_dir/norm(noise_dir);
                    uv_with_noise = uv+normrnd(0, noise_level)*noise_dir;
                    if (abs(uv_with_noise(1)) > img_width/2 || abs(uv_with_noise(2)) > img_height/2)
                        continue;
                    end
                    break;
                end
            else
                while(true)
                    uv_with_noise = [img_width*(rand(1)-0.5); img_height*(rand(1)-0.5)];
                    noise_magnitude = norm(uv-uv_with_noise);
                    if (noise_magnitude < outlier_min_noise)
                        continue;
                    end
                    break;
                end
            end

            break;
        end

        R{i} = R_temp;
        t{i} = t_temp;
        uv_measured{i} = uv_with_noise;
    end
end


function [cam_pos, P, A, M3, M4, M5] = Precompute(n_cameras, K, R, t)


    cam_pos = cell(1,n_cameras);
    P = cell(1,n_cameras);
    A = cell(1,n_cameras); % necessary for computing the Jacobian
    M3 = nan(4,n_cameras);
    M4 = nan(4,n_cameras);
    M5 = nan(4,n_cameras);
    
    for i = 1:n_cameras
        r11 = R{i}(1,1); r12 = R{i}(1,2); r13 = R{i}(1,3);
        r21 = R{i}(2,1); r22 = R{i}(2,2); r23 = R{i}(2,3);
        r31 = R{i}(3,1); r32 = R{i}(3,2); r33 = R{i}(3,3);
        tx = t{i}(1); ty = t{i}(2); tz = t{i}(3);
        k11 = K{i}(1,1); k12 = K{i}(1,2); k21 = K{i}(2,1); k22 = K{i}(2,2);

        b1 = r11*[0;r32;r33;tz]-r31*[0;r12;r13;tx];
        b2 = r21*[0;r32;r33;tz]-r31*[0;r22;r23;ty];
        b3 = r12*[r31;0;r33;tz]-r32*[r11;0;r13;tx];
        b4 = r22*[r31;0;r33;tz]-r32*[r21;0;r23;ty];
        b5 = r13*[r31;r32;0;tz]-r33*[r11;r12;0;tx];
        b6 = r23*[r31;r32;0;tz]-r33*[r21;r22;0;ty];
        a1 = k11*b1+k12*b2;
        a2 = k21*b1+k22*b2;
        a3 = k11*b3+k12*b4;
        a4 = k21*b3+k22*b4;
        a5 = k11*b5+k12*b6;
        a6 = k21*b5+k22*b6;
        A{i} = [a1,a2,a3,a4,a5,a6];
        
        cam_pos{i} = -R{i}'*t{i};
        P{i} = [R{i}, t{i}];
        
        M3(:,i) = K{i}(1,1)*P{i}(1,:)'+K{i}(1,2)*P{i}(2,:)';
        M4(:,i) = K{i}(2,1)*P{i}(1,:)'+K{i}(2,2)*P{i}(2,:)';
        M5(:,i) = P{i}(3,:)';
    end

end



function [output_pairs] = GetSamplePairs(input_vec, max_n_samples, bool_random_forced)
    n_sample = length(input_vec);
    possible_n_pairs = n_sample*(n_sample-1)/2;
    if(~bool_random_forced && possible_n_pairs <= max_n_samples)
        output_pairs = zeros(2,possible_n_pairs);
        c = 0;
        for i = input_vec
            for j = input_vec
                if (i<=j)
                    continue;
                end
                c = c + 1;
                output_pairs(:,c) = [i,j];
            end
        end
    else
        max_n_samples = min(max_n_samples, possible_n_pairs);
        output_pairs = zeros(2,max_n_samples);
        idx_shuffled = input_vec(randperm(n_sample));
        c = 0;
        gap = 1;
        sample1_idx = 0;
        while (c<max_n_samples)
            sample1_idx = sample1_idx + 1;
            sample2_idx = sample1_idx+gap;
            if (sample1_idx > n_sample || sample2_idx > n_sample)
                sample1_idx = 0;
                gap = gap + 1;
                continue;
            end
            sample1 = idx_shuffled(sample1_idx);
            sample2 = idx_shuffled(sample2_idx);
            c = c + 1;
            output_pairs(:,c) = [sample1,sample2];
        end
    end
end

function error_2D = GetMean2DErrorOfTrueInliers(x, n_cameras, n_inlier, K, P, uv_measured)

    M1 = nan(1,n_cameras);
    M2 = nan(1,n_cameras);
    M3 = nan(4,n_cameras);
    M4 = nan(4,n_cameras);
    M5 = nan(4,n_cameras);
    for i = 1:n_cameras
        M1(1,i) = K{i}(1,3)-uv_measured{i}(1);
        M2(1,i) = K{i}(2,3)-uv_measured{i}(2);
        M3(:,i) = K{i}(1,1)*P{i}(1,:)'+K{i}(1,2)*P{i}(2,:)';
        M4(:,i) = K{i}(2,1)*P{i}(1,:)'+K{i}(2,2)*P{i}(2,:)';
        M5(:,i) = P{i}(3,:)';
    end
    
    M6 = [x;1]'*M5;
    M7 = M1+([x;1]'*M3)./M6;
    M8 = M2+([x;1]'*M4)./M6;
    squared_reproj_error_all = M7.*M7+M8.*M8;
    
    error_2D = 0;
    c = 0;
    for i = 1:n_inlier
        if (M6(i) <= 0)
            continue;
        end
        error_2D = error_2D + sqrt(squared_reproj_error_all(i));
        c = c + 1;
    end
    error_2D = error_2D/c;
    
    if (error_2D == 0)
        error_2D = nan;
    end
end

function error_2D = GetMean2DErrorOfEstimatedInliers(x, M1, M2, M3, M4, M5, idx_inlier)
    M6 = [x;1]'*M5;
    M7 = M1+([x;1]'*M3)./M6;
    M8 = M2+([x;1]'*M4)./M6;
    squared_reproj_error_all = M7.*M7+M8.*M8;
    
    error_2D = 0;
    c = 0;
    for i = idx_inlier
        if (M6(i) <= 0)
            continue;
        end
        error_2D = error_2D + sqrt(squared_reproj_error_all(i));
        c = c + 1;
    end
    error_2D = error_2D/c;
    
    if (error_2D == 0)
        error_2D = nan;
    end
end


function [min_level_of_degeneracy] = GetMinLevelOfDegeneracy(cam_pos, ray_w, K, uv_measured, R, idx, n_sample_mLoD)
    
    min_level_of_degeneracy = inf;
    idx = reshape(idx, [1, numel(idx)]);
    
    sample_pairs = GetSamplePairs(idx, n_sample_mLoD, false);

    for k = 1:size(sample_pairs,2)
        i = sample_pairs(1,k);
        j = sample_pairs(2,k);

        ray_i = ray_w(:,i);
        if (ray_i(1)==0 && ray_i(2)==0 && ray_i(3)==0)
            ray_i = K{i}\[uv_measured{i};1];
            ray_i = ray_i/norm(ray_i);
            ray_i = R{i}'*ray_i;
            ray_w(:,i) = ray_i;
        end

        ray_j = ray_w(:,j);
        if (ray_j(1)==0 && ray_j(2)==0 && ray_j(3)==0)
            ray_j = K{j}\[uv_measured{j};1];
            ray_j = ray_j/norm(ray_j);
            ray_j = R{j}'*ray_j;
            ray_w(:,j) = ray_j;
        end
    
        t_ij = cam_pos{i}-cam_pos{j};
        t_ij = t_ij/norm(t_ij);
        f_i = ray_i;
        f_j = ray_j;
        
        level_of_degeneracy = max(abs([t_ij'*f_i, t_ij'*f_j, f_i'*f_j]));

        if (level_of_degeneracy < min_level_of_degeneracy)
            min_level_of_degeneracy = level_of_degeneracy;
        end
    end
  
end



function [x_est, idx_consensus, max_n_consensus, min_cost, k_max, ray_w] = TwoViewRansac(j, k, x_est, idx_consensus, max_n_consensus, min_cost, k_max, confidence, n_cameras, epipolar_thr, cos_angle_thr, reproj_thr, M1, M2, M3, M4, M5, K, R, P, cam_pos, uv_measured,ray_w)
    
    ray_j = ray_w(:,j);
    if (ray_j(1)==0 && ray_j(2)==0 && ray_j(3)==0)
        ray_j = K{j}\[uv_measured{j};1];
        ray_j = ray_j/norm(ray_j);
        ray_j = R{j}'*ray_j;
        
        ray_w(:,j) = ray_j;
    end
    
    ray_k = ray_w(:,k);
    if (ray_k(1)==0 && ray_k(2)==0 && ray_k(3)==0)
        ray_k = K{k}\[uv_measured{k};1];
        ray_k = ray_k/norm(ray_k);
        ray_k = R{k}'*ray_k;
        
        ray_w(:,k) = ray_k;
    end

        
        
 
    % Epipolar error should be below threshold
    t_jk = cam_pos{j}-cam_pos{k};
    t_norm = norm(t_jk);
    t_jk = t_jk/t_norm;
    cross_fj_fk = [...
        ray_j(2)*ray_k(3)-ray_j(3)*ray_k(2);...
        ray_j(3)*ray_k(1)-ray_j(1)*ray_k(3);...
        ray_j(1)*ray_k(2)-ray_j(2)*ray_k(1)];
    e_epipolar = abs(t_jk'*cross_fj_fk);   

    if (e_epipolar > epipolar_thr)
        return;
    end

    % Parallax must be below 90 deg, above threshold. 
    p_jk = ray_j'*ray_k;
    if (p_jk < 0 || p_jk > cos_angle_thr)
        %disp(['Skip I:' num2str(i), ', ', num2str(j)])
        return;
    end

    % Point must be far from the epipoles
    q_jk = t_jk'*ray_j;
    r_jk = t_jk'*ray_k;
    if (abs(q_jk) > cos_angle_thr || abs(r_jk) > cos_angle_thr)
        %disp(['Skip II:' num2str(i), ', ', num2str(j)])
        return;
    end

    % Midpoint should satisfy the cheirality.
    tau_i = p_jk*r_jk-q_jk;
    tau_j = -p_jk*q_jk+r_jk;
    if (tau_i < 0 || tau_j < 0)
        %disp(['Skip III:' num2str(i), ', ', num2str(j)])
        return;
    end

    % Obtain the midpoint.
    s_jk = t_norm/(1-p_jk'*p_jk);
    lambda_j = s_jk*tau_i;
    lambda_k = s_jk*tau_j;
    x_anchor_j = cam_pos{j}+lambda_j*ray_j;
    x_anchor_k= cam_pos{k}+lambda_k*ray_k;
    x_mid = 0.5*(x_anchor_j+x_anchor_k);

    
    % Cheirality
    
    x_j = P{j}*[x_mid;1];
    x_k = P{k}*[x_mid;1];

    if (x_j(3) < 0 || x_k(3) < 0)
        %disp(['Skip III:' num2str(i), ', ', num2str(j)])
        return;
    end
    
    % Small reprojection error
    u_j = uv_measured{j}(1); v_j = uv_measured{j}(2);
    u_k = uv_measured{k}(1); v_k = uv_measured{k}(2);
    
    uv_reproj_j = K{j}*x_j/x_j(3);
    uv_error_j = [u_j;v_j;1]-uv_reproj_j;
    reproj_error_j = uv_error_j'*uv_error_j;
    if (reproj_error_j > reproj_thr)
        return;
    end

    uv_reproj_k = K{k}*x_k/x_k(3);
    uv_error_k = [u_k;v_k;1]-uv_reproj_k;
    reproj_error_k = uv_error_k'*uv_error_k;
    if (reproj_error_k > reproj_thr)
        return;
    end
    
    % Count the consensus.
    M6 = [x_mid;1]'*M5;
    M7 = M1 + ([x_mid;1]'*M3)./M6;
    M8 = M2 + ([x_mid;1]'*M4)./M6;
    squared_reproj_error_all = M7.*M7+M8.*M8;
    squared_reproj_error_all(M6<=0) = inf;

    cost = 0;
    for ii = 1:length(squared_reproj_error_all)
        cost = cost + min(squared_reproj_error_all(ii), reproj_thr);
    end

    bool_consensus = squared_reproj_error_all<reproj_thr;

    if (cost < min_cost)
        min_cost = cost;

        bool_max_consensus = bool_consensus;
        max_n_consensus = sum(bool_consensus);

        x_est = x_mid;
        idx_consensus = find(bool_max_consensus);

        inlier_ratio_est = (max_n_consensus+1)/n_cameras;
        k_max = log(1-confidence)/log(1-inlier_ratio_est^2);
        %disp([num2str(k), ', ', num2str(k_max)])

    end
end

function [x_est, idx_consensus, max_n_consensus] = DLT(x_est, idx_consensus, n_cameras, reproj_thr, M1, M2, M3, M4, M5, K, P, uv_measured)

    
    bool_consensus = zeros(1,n_cameras);
    
    for it = 1:10
        M6 = [x_est;1]'*M5;
        M7 = M1+([x_est;1]'*M3)./M6;
        M8 = M2+([x_est;1]'*M4)./M6;
        squared_reproj_error_all = M7.*M7+M8.*M8;
        squared_reproj_error_all(M6<=0) = inf;

        
        bool_consensus_prev = bool_consensus;
        bool_consensus = squared_reproj_error_all<reproj_thr;
        
        if (sum(abs(bool_consensus_prev-bool_consensus))==0 || sum(bool_consensus)<2)
            break;
        end
        
        A = [];
        for i = 1:n_cameras
            if (bool_consensus(i))
                Q = K{i}*P{i};
                u = uv_measured{i}(1);
                v = uv_measured{i}(2);
                row1 = u*Q(3,:)-Q(1,:);
                row2 = v*Q(3,:)-Q(2,:);

                A = [A; row1; row2];
            end
        end
        [~,S_,V_] = svd(A);
        [~, min_idx] = min(diag(S_));
        x_est = V_(:, min_idx);
        x_est = x_est/x_est(end);
        x_est = x_est(1:3);
    end
    
    % Count the consensus.
    M6 = [x_est;1]'*M5;
    M7 = M1+([x_est;1]'*M3)./M6;
    M8 = M2+([x_est;1]'*M4)./M6;
    squared_reproj_error_all = M7.*M7+M8.*M8;
    squared_reproj_error_all(M6<=0) = inf;

    bool_consensus = squared_reproj_error_all<reproj_thr;
    idx_consensus = find(bool_consensus);
    max_n_consensus = sum(bool_consensus);
end

function [x_est, idx_consensus, max_n_consensus] = GaussNewton(x_est, idx_consensus, A, P, n_cameras, reproj_thr, M1, M2, M3, M4, M5)

    bool_consensus = zeros(1,n_cameras);
    mean_reproj_error = 0;
    c = 0;
    
    x_gn = x_est;
    for it = 1:10

        M6 = [x_gn;1]'*M5;
        M7 = M1+([x_gn;1]'*M3)./M6;
        M8 = M2+([x_gn;1]'*M4)./M6;
        residuals = [M7;M8];
        squared_reproj_error_all = M7.*M7+M8.*M8;
        squared_reproj_error_all(M6<=0) = inf;
        
        mean_reproj_error_prev = mean_reproj_error;
        mean_reproj_error = mean(sqrt(squared_reproj_error_all(squared_reproj_error_all<reproj_thr)));
        
        bool_consensus_prev = bool_consensus;
        bool_consensus = squared_reproj_error_all<reproj_thr;
        bool_consensus_changed = sum(abs(bool_consensus-bool_consensus_prev)) > 0;
        
        if (bool_consensus_changed)
            c = 0;
        else
            c = c +1;
        end
        
        if (sum(bool_consensus) < 2)
            break;
        end
        
        if (c == 2 && abs(mean_reproj_error-mean_reproj_error_prev) < 0.1)
            break;
        end
        
        residuals_new = [];
        J_est = [];

        for i = 1:n_cameras
            if (bool_consensus(i))
                residuals_new = [residuals_new; residuals(:,i)];

                J_temp = A{i}'*[x_gn;1];
                J_temp = reshape(J_temp, [2,3])/((P{i}(3,:)*[x_gn;1])^2); 
                J_est = [J_est; J_temp];
            end
        end
        
        update = -pinv(J_est)*residuals_new;
        
        x_gn = x_gn+update;
        
    end

    x_est = x_gn;
    idx_consensus = find(bool_consensus);
    max_n_consensus = sum(bool_consensus);
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





function [x_est, uncertainty_3D] = RunMultiviewTriangulation(method, uncertainty, cam_pos, R,P,K,uv_measured,A, M3, M4, M5)

    %% 1. Initialize
    
    uncertainty_3D = nan;
    x_est = nan;
    idx_consensus = [];
    
    n_cameras = numel(cam_pos);
    
    % Construct matrix M1 and M2 
    M1 = nan(1,n_cameras);
    M2 = nan(1,n_cameras);
    for i = 1:n_cameras
        M1(1,i) = K{i}(1,3)-uv_measured{i}(1);
        M2(1,i) = K{i}(2,3)-uv_measured{i}(2);
    end


    ray_w = zeros(3,n_cameras);

    
    
    %% 2. Run two-view RANSAC

    confidence = 0.99;
    cos_angle_thr = cosd(4);
    reproj_thr = 10^2;
    epipolar_thr = 0.01;
    
    n_sample = min(10000, n_cameras*(n_cameras-1)/2);
    n_sample_mLoD = 100;
    n_consensus_thr = 2;
    
    min_cost = inf;
    inlier_ratio_est = 3/n_cameras;
    k_max = log(1-confidence)/log(1-inlier_ratio_est^2);
    max_n_consensus = 0;

    sample_pairs = GetSamplePairs(1:n_cameras, inf, true);

    for k = 1:size(sample_pairs,2)
        if (k == n_sample)
            break;
        end
        if (k > k_max && max_n_consensus >= n_consensus_thr)
            break;
        end
        i = sample_pairs(1,k);
        j = sample_pairs(2,k);
        [x_est, idx_consensus, max_n_consensus, min_cost, k_max, ray_w] = TwoViewRansac(i, j, x_est, idx_consensus, max_n_consensus, min_cost, k_max, confidence, n_cameras, epipolar_thr, cos_angle_thr, reproj_thr, M1, M2, M3, M4, M5, K, R, P, cam_pos, uv_measured,ray_w);

    end

    if (max_n_consensus < n_consensus_thr)
        disp('Not enough inliers!')
        return;
    end
    
    
    %% 3. Refine the initial solution and the inlier set:
    switch method
        case 'GN'
            [x_est, idx_consensus, max_n_consensus] = GaussNewton(x_est, idx_consensus, A, P, n_cameras, reproj_thr, M1, M2, M3, M4, M5);
        case 'DLT'
            [x_est, idx_consensus, max_n_consensus] = DLT(x_est, idx_consensus, n_cameras, reproj_thr, M1, M2, M3, M4, M5, K, P, uv_measured);
    end
    
    %% 4. Estimate the 3D uncertainty:
    error_2D_estimated_inliers = GetMean2DErrorOfEstimatedInliers(x_est, M1, M2, M3, M4, M5, idx_consensus);
    mLoD_approx = GetMinLevelOfDegeneracy(cam_pos, ray_w, K, uv_measured, R, idx_consensus, n_sample_mLoD);
    uncertainty_3D = TrilinearInterpolation(acosd(mLoD_approx), error_2D_estimated_inliers, max_n_consensus, uncertainty);
                        

end

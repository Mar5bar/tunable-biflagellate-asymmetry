function h = draw_robot(position, phi, t, ind, params)
    hold on
    h = []; 
    position = [position(:);0];
    orientation = rot_mat(phi);
    % Draw the body.
    bodyX = params.r*cos(linspace(0,2*pi,1e2));
    bodyY = params.r*sin(linspace(0,2*pi,1e2));
    for i = 1 : numel(bodyX)
        mappedVec = orientation * [bodyX(i); bodyY(i); 0] + position;
        bodyX(i) = mappedVec(1);
        bodyY(i) = mappedVec(2);
    end
    h = [h; patch(bodyX, bodyY, 0.7*[1,1,1])];

    % Draw the legs.
    alphas = params.alphas(ind,:);
    betas = params.betas(ind,:);
    for i = 1 : 2
        n = compute_n(i,params);
        t = compute_t(i,params);
        alpha = alphas(i);
        beta = betas(i);
        % Compute everything in the swimmer frame.
        point = attachmentPoint(i, params);
        pPar = cos(alpha)*n + sin(alpha)*t;
        qPar = cos(alpha+beta)*n + sin(alpha+beta)*t;
        midpoint = point + params.l1 * pPar;
        endpoint = midpoint + params.l2 * qPar;
        % Map to the lab frame.
        point = orientation * point + position;
        midpoint = orientation * midpoint + position;
        endpoint = orientation * endpoint + position;
        coords = [point, midpoint, endpoint];

        h = [h; plot(coords(1,:),coords(2,:),'black')];
    end
    axis equal
end
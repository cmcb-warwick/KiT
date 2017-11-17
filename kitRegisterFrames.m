function dataStruct=kitRegisterFrames(job,reader,dataStruct)
%KITREGISTERFRAMES Rotates spot coordinates using image registration

md = job.metadata;
options = job.options;

planeFit(1:md.nFrames) = struct('plane',[],'planeCoord',[],'planeVectorClassifier', 0, ...
    'planeVectors',[],'planeOrigin',[],'angle',[],'translation',[],...
    'inlierIdx',[],'unalignedIdx',[],'laggingIdx',[],'phase','e',...
    'tform',[]);
initCoord = dataStruct.initCoord;

% Setup image registration optimizer.
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 300;
optimizer.MinimumStepLength = 5e-4;

% Read frame by frame.
maxMergeChannels = 3;
for t=1:md.nFrames
    cImg = ones(job.cropSize([2 1]));
    for c=1:min(md.nChannels, maxMergeChannels)
        % Read stack.
        img = kitReadImageStack(reader, md, t, c, job.crop, 0);
        % Max project.
        img = max(img, [], 3);
        % Flip image.
        img = img';
        % Enhance contrast.
        img = imadjust(img, stretchlim(img,0),[]);
        % Multiply blend.
        cImg = immultiply(cImg,img);
    end

    if t>1
        % Derive rigid transformation
        tform = imregtform(cImg,firstImg,'rigid',optimizer,metric);
        % Transform image for next frame.
        transImg = imwarp(cImg,(tform));
        % Extract affine transform matrix.
        tMat = tform.T;
        % Convert to microns.
        tMat(3,1:2) = tMat(3,1:2) .* md.pixelSize(1:2);
        tform.T = tMat;
        planeFit(t).tform = tform;
    else
        tMat = eye(3);
        planeFit(t).tform = affine3d();
        transImg = cImg;
        firstImg = cImg;
    end
    planeFit(t).angle = asin(1-tMat(1,1))*180/pi;
    planeFit(t).translation = [tMat(3,1:2) 0];

    nCoords = size(initCoord(t).allCoord, 1);
    % Transform spots as homogeneous coordinates.
    homCoords = [initCoord(t).allCoord(:,1:2) ones(nCoords,1)];
    planeFit(t).rotatedCoord = homCoords*tMat;
    % Error propagation.
    errorPropMat = tMat.^2;
    homCoords = [initCoord(t).allCoord(:,4:5) ones(nCoords,1)];
    planeFit(t).rotatedCoord(:,4:6) = sqrt(homCoords.^2*errorPropMat);
    % Replace homogenous coordinate with z.
    planeFit(t).rotatedCoord(:,[3 6]) = initCoord(t).allCoord(:,[3 6]);
    planeFit(t).planeOrigin = planeFit(t).translation;
    planeFit(t).planeVectors = [[tMat(1:2,1:2) [0; 0]]; 0 0 1];

    if options.debug.showPlaneFit
        figure(1);
        imshowpair(firstImg,transImg);
        hold on;
        po = planeFit(t).planeOrigin(1:2) ./ md.pixelSize(1:2);
        px = tMat(1,1:2);
        py = tMat(2,1:2);
        plength = 40;
        p1 = [po; po+plength*px];
        p2 = [po; po+plength*py];
        plot(po(1),po(2),'wx');
        plot(p1(:,1), p1(:,2), '-r'); % X
        plot(p2(:,1), p2(:,2), '-c'); % Y
        hold off;
    end
    
    prevImg = transImg;
end


dataStruct.planeFit = planeFit;
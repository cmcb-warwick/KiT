function kitMakeZxProj(job,plotPlane)
%KITMAKEZXPROJ Renders a movie in Z-X projection

if nargin<3
    useRotated = 0;
end

vFilename = strrep(job.output, '.mat', '_zxproj.mj2');

% Open movie.
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie));

vWriter = VideoWriter(vFilename, 'Motion JPEG 2000');
vWriter.FrameRate = 5;
open(vWriter);

% Read frame by frame.
levels = 2^(md.nBytesPerPixel*8);
maps = {@green, @red, @blue};
markers = 'xo+*sd';
maxMergeChannels = 3;
for i=1:md.nFrames
    rgbImg = zeros([job.cropSize([1 3]), 3]);
    %rgbImg = zeros([job.cropSize([2 1]), 3]);
    for c=1:min(md.nChannels, maxMergeChannels)
        % Read stack.
        img = kitReadImageStack(reader, md, i, c, job.crop, 0);

        % Max project.
        img = squeeze(max(img, [], 2));
        % Flip image.
        %img = img';

        % Merge into RGB image.
        rgbImg = imadd(ind2rgb(imadjust(img, stretchlim(img,0),[]), ...
                     maps{c}(levels)), rgbImg);
    end


    figure(1);
    imshow(rgbImg);
    axis square;
    %set(gca,'ydir','normal');

    c=job.options.coordSystemChannel;
    % Plot plane fit.
    if plotPlane && ~isempty(job.dataStruct{c}.planeFit(i).planeVectors)
        px = [job.dataStruct{c}.planeFit(i).planeVectors([1 3],1)' 0];
        pz = [job.dataStruct{c}.planeFit(i).planeVectors([1 3],3)' 0];
        po = [job.dataStruct{c}.planeFit(i).planeOrigin([1 3])./job.metadata.pixelSize([1 3]) 1];

        hold on
        plength = 40;
        pzlength = 1;
        p1 = [po; po+plength*px];
        p2 = [po; po+pzlength*pz];
        plot(p1(:,2), p1(:,1), '-r');
        plot(p2(:,2), p2(:,1), '-c');
        plot(po(2), po(1), 'yx');

        hold off
    end

    % Save frame.
    writeVideo(vWriter, getframe);
end
close(vWriter);

%% NESTED FUNCTIONS
function rgb=green(n)
    x=linspace(0,1,n);
    rgb=[zeros(n,1), x', zeros(n,1)];
end

function rgb=blue(n)
    x=linspace(0,1,n);
    rgb=[zeros(n,1), zeros(n,1), x'];
end

function rgb=red(n)
    x=linspace(0,1,n);
    rgb=[x', zeros(n,1), zeros(n,1)];
end

end
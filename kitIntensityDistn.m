function m=kitIntensityDistn(job,reader,metadata,channel,z,frames,dynamic,crop,nomean)

if nargin<5
    z = [];
end
if nargin<6
    frames = [];
end
if nargin<8
    crop = [];
end
if nargin<9
    nomean = 0;
end


% Open movie.
if isstruct(job)
    movieFile = job.movie;
elseif ischar(job)
    movieFile = job;
else
    [filename,path] = uigetfile([],'Select movie file');
    movieFile = fullfile(path,filename);
end

if isempty(frames)
    frames=1:metadata.nFrames;
end

if dynamic
    if nomean && isempty(crop)
        m=zeros(prod(job.metadata.frameSize),length(frames));
    else
        m=zeros(length(frames),1);
    end
else
    m=[];
end

for t = frames
    if isempty(z)
        stack = kitReadImageStack(reader, metadata, t, channel, crop, 0);
    else
        stack = kitReadImagePlane(reader, metadata, t, channel, z, crop, 0);
    end
    if dynamic
        if nomean && isempty(crop)
            m(:,t) = stack(:);
        else
            m(t) = mean(stack(:));
        end
    else
        m = [m; stack(:)];
    end
end

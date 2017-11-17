function filters=createFilters(nDims,dataProperties)

backFilterParms = dataProperties.FT_SIGMA * 14; % 14 is 15-1
backFilterParms(4:6) = roundOddOrEven(backFilterParms,'odd','inf');
% Constrain filter by number of planes in stack.
backFilterParms(6) = min(backFilterParms(6), (dataProperties.movieSize(3)-1)*2+1);
signalFilterParms = dataProperties.FILTERPRM;

switch nDims
  case 3
    backgroundFilter = GaussMask3D(backFilterParms(1:3),...
        backFilterParms(4:6),[],1,[],[],1);
    signalFilter = GaussMask3D(signalFilterParms(1:3),...
        signalFilterParms(4:6),[],1,[],[],1);
    borderMode = 2;
    % make separated noise mask
    noiseMask = {...
        ones(signalFilterParms(4),1,1)./signalFilterParms(4),...
        ones(1,signalFilterParms(5),1)./signalFilterParms(5),...
        ones(1,1,signalFilterParms(6))./signalFilterParms(6),...
        };
  case 2
    backgroundFilter = GaussMask2D(backFilterParms(1:2),...
        backFilterParms(4:5),[],1,[]);
    signalFilter = GaussMask2D(signalFilterParms(1:2),...
        signalFilterParms(4:5),[],1,[]);
    borderMode = 1;
    % make separated noise mask
    noiseMask = {...
        ones(signalFilterParms(4),1,1)./signalFilterParms(4),...
        ones(1,signalFilterParms(5),1)./signalFilterParms(5)};
  otherwise
    error('Unsupported dimensions: %d',ndims);
end

filters.background = backgroundFilter;
filters.signal = signalFilter;
filters.backgroundP = backFilterParms;
filters.signalP = signalFilterParms;
filters.border = borderMode;
filters.noise = noiseMask;

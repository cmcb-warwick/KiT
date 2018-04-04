function spotSel = combineSpotSelections(spotSels)
% Use like newSpotSel = combineSpotSelections({spotSel_1,spotSel_2,...});
%
% Copyright (c) 2017 C. A. Smith

if ~iscell(spotSels)
    error('Please provide more than one spot selection structure in cell format.')
end
nExpts = length(spotSels);

spotSel = spotSels{1};
method = spotSel.dataType;
kitLog('Combining %i spot selection structures made using %s.',nExpts,method)

for iExpt = 2:nExpts
    if ~strcmp(spotSels{iExpt}.dataType,method)
        kitLog('Experiment %i ignored: contains data selected using %s.',iExpt,method);
        continue
    else
        spotSel.selection{iExpt} = spotSels{iExpt}.selection{1};
        if isfield(spotSels{iExpt},'rawSelection')
            spotSel.rawSelection{iExpt} = spotSels{iExpt}.rawSelection{1};
        end
    end
end

end
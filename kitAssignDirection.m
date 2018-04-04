function [direction,vote]=kitAssignDirection(normalCoord,options)

% get specific options
mode = options.assignMode;
plot=0;
exp = options.assignExpWeight;
cohSteps = options.minConsSteps;
buffer = options.switchBuffer;

% Standardize direction: +ve means P, -ve mean AP.
% Compute displacement.
displ = diff(normalCoord);

% Designate P as 1 and AP as -1.
fwdback = (displ > 0) - (displ < 0);

% Set steps not bracketed by same step to indeterminate 0.
direction = zeros(size(fwdback));
n=length(direction);

switch mode
  case 'voting'
    
    w=3;
    if exp
      % Exponentially weight displacements.
      lambda = expfit(abs(displ(~isnan(displ))));
      displ = expcdf(abs(displ),lambda);
      fwdback = displ.*fwdback;
      vcut=lambda;
    else
      vcut=0.25;
    end
    
    vote=zeros(size(direction));
    for i=1:n
        range=max(1,i-w):min(n,i+w);

        if exp
          v=nanmean(fwdback(range));
        else
          v=mean(fwdback(range));
        end
        if v>=vcut
            direction(i) = 1;
        elseif v<=-vcut
            direction(i) = -1;
        end
        vote(i)=v;
    end

    if plot
      figure;
      t=1:length(normalCoord);
      m=nanmean(normalCoord);
      plot(t,normalCoord-m,'b+-',t(direction==1),normalCoord(direction==1)-m,'rx',...
           t(direction==-1),normalCoord(direction==-1)-m,'go',t,[vote; nan],'r-');
    end

  case 'absolute'
    
    % loop over coordinates
    for iCoord = 1:n-cohSteps+1

      % if all steps in a given direction, set as that direction
      if abs(sum(fwdback(iCoord:iCoord+cohSteps-1))) == cohSteps
        direction(iCoord+buffer:iCoord+cohSteps-1-buffer) = fwdback(iCoord);
      end
      
    end
  
  otherwise
    error('Unknown direction assignment mode');
end


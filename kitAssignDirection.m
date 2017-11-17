function [direction,vote]=kitAssignDirection(normalCoord,varargin)

opts.mode='voting';
opts.plot=0;
opts.exp=0;
opts=processOptions(opts,varargin{:});

switch opts.mode
  case 'voting'
    % Standardize direction: +ve means P, -ve mean AP.
    % Compute displacement.
    displ = diff(normalCoord);

    % Designate P as 1 and AP as -1.
    fwdback = (displ > 0) - (displ < 0);

    % Set steps not bracketed by same step to indeterminate 0.
    direction = zeros(size(fwdback));
    w=3;
    if opts.exp
      % Exponentially weight displacements.
      lambda = expfit(abs(displ(~isnan(displ))));
      displ = expcdf(abs(displ),lambda);
      fwdback = displ.*fwdback;
      vcut=lambda;
    else
      vcut=0.25;
    end
    n=length(direction);
    vote=zeros(size(direction));
    for i=1:n
        range=max(1,i-w):min(n,i+w);

        if opts.exp
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

    if opts.plot
      figure;
      t=1:length(normalCoord);
      m=nanmean(normalCoord);
      plot(t,normalCoord-m,'b+-',t(direction==1),normalCoord(direction==1)-m,'rx',...
           t(direction==-1),normalCoord(direction==-1)-m,'go',t,[vote; nan],'r-');
    end

  otherwise
    error('Unknown direction assignment mode');
end


function h=label(p,s,pos,sz)
% Create text label
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','text','String',s,'Units','characters','Position',pos,'FontSize',sz,'HorizontalAlignment','left');
end

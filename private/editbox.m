function h=editbox(p,s,pos,sz)
% Create edit box
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','edit','String',s,'Units','characters','Position',pos,'FontSize',sz,'HorizontalAlignment','left');
end

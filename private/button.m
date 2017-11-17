function h=button(p,s,pos,cb,sz)
% Create button
  if nargin<4
    cb = '';
  end
  if nargin<5
    sz=12;
  end
  h = uicontrol(p,'String',s,'Units','characters','Position',pos,'FontSize',sz,'Callback',cb);
end

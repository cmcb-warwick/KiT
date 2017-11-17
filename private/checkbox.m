function h=checkbox(p,s,pos,cb,sz)
% Create checkbox
  if nargin<4
    cb = '';
  end
  if nargin<5
    sz=12;
  end
  h = uicontrol(p,'Style','checkbox','String',s,'Units','characters','Position',pos,'FontSize',sz,'Callback',cb);
end

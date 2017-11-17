function h=popup(p,s,pos,cb,sz)
% Create popup menu
  if nargin<4
    cb = '';
  end
  if nargin<5
    sz=12;
  end
  h = uicontrol(p,'Style','popupmenu','String',s,'Units','characters','Position',pos,'FontSize',sz,'Callback',cb);
end

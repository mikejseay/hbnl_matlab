function [] = GUI_plotparams(pp)

%get the input structure's field names
field_name=fieldnames(pp);

%ignore certain fields
fields2ignore={'chosen_freq','chosen_cond,'chosen_s','
fieldnames( strcmp(fieldnames,'

n_fields=length(field_name);
height=n_fields*30;

%create the overall figure
S.fh = figure('units','pixels',...
              'position',[300 300 600 height],...
              'menubar','none',...
              'numbertitle','off',...
              'name','Plotting Parameters',...
              'resize','off');


%for each field
for field=1:n_fields
    %label to the left
S.lbl{field} = uicontrol('style','text',...
                 'units','pixels',...
                 'position',[10 50 280 30],...
                 'fontsize',14,...
                 'string',field_name{field});
     %interpret the structure contents
     field_contents=getfield(pp,field_name);
     if iscell(field_contents)
     elseif ischar(field_contents)
     elseif 
S.ed{field} = uicontrol('style','edit',...
                 'units','pixels',...
                 'position',[300 50 280 30],...
                 'fontsize',14,...
                 'string',field_contents);
end

%put in the set params button
S.pb = uicontrol('style','push',...
                 'units','pixels',...
                 'position',[240 10 100 30],...
                 'fonts',14,...
                 'str','Set Params',...
                 'callback',{@pb_call,S});

             
function [] = pb_call(varargin)
% Callback for the button labeled PUSH_1.
S = varargin{3}; % Get the structure.
AX = findobj('type','axes');  % Look for existing axes objects.

if isempty(AX)
   % The GUI could popup a message box saying, "No axes to title!" or 
   % something similar, then return from here.
   figure
   AX = axes;
end

title(AX(1),get(S.ed,'string')) % The first one in the list was current!


    









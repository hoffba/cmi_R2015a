function killWaitbars
set(0,'ShowHiddenHandles','on');
h = get(0,'Children');
ht = get(h,'Tag');
i = strcmp('TMWWaitbar',ht);
delete(h(i));
set(0,'ShowHiddenHandles','off');
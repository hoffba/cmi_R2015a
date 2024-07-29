function plot_graph(ha,n,l)
% given structures n of nodes and l of links from Skel2Graph3D, plot links
hold(ha,'on');
s=size(l);
for i=1:s(2)
    n1=l(i).n1; n2=l(i).n2;
    n1=[n(n1).comx n(n1).comy n(n1).comz];
    n2=[n(n2).comx n(n2).comy n(n2).comz];
    plot3(ha,[n1(1) n2(1)],[n1(2) n2(2)],[n1(3) n2(3)],'black-');
    plot3(ha,n1(1),n1(2),n1(3),'ro');
    plot3(ha,n2(1),n2(2),n2(3),'ro');
end
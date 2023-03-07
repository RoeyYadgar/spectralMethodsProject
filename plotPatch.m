function plotPatch(patch,pos)
%plots a patch
%inputs:
%   patch - patch graph 
%   pos - position of nodes of the patch

%if pos was not given as an input, use the Pos attribute of the graph
if(nargin == 1)
   pos = patch.Nodes.Pos; 
end

%if pos is Nx2 or Nx3 matrix use each column as the cooridnates of x,y,z
if(size(pos,2) >= 2)
    xPos = pos(:,1);
    yPos = pos(:,2);
    if(size(pos,2) == 3)
        zPos = pos(:,3);
    else
        zPos = zeros(size(xPos));
    end
end
%if pos is complex use the real and imaginary parts as the cooridnates of x,y
if(~isreal(pos))
    xPos = real(pos);
    yPos = imag(pos);
    zPos = zeros(size(xPos));

end

%plot the graph
p = plot(patch,'XData',xPos,'YData',yPos,'ZData',zPos);
p.LineWidth = 2;
view(0,90)

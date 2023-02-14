function plotPatch(patch,pos)


if(nargin == 1)
   pos = patch.Nodes.Pos; 
end
if(size(pos,2) >= 2)
    xPos = pos(:,1);
    yPos = pos(:,2);
    if(size(pos,2) == 3)
        zPos = pos(:,3);
    else
        zPos = zeros(size(xPos));
    end
end
if(~isreal(pos))
    xPos = real(pos);
    yPos = imag(pos);
    zPos = zeros(size(xPos));

end

    
plot(patch,'XData',xPos,'YData',yPos,'ZData',zPos)
view(0,90)

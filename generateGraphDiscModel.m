function G = generateGraphDiscModel(points,rho,eta)
%Generate a graph using the Disc model assumption, i.e. (between any pair
%of points there is an edge iff their distance is smaller then rho) 
%inputs:
%   points - N x 2 matrix containing the coordiantes of the points
%   rho - distance threshold of the disc model
%   eta - noise level
%outputs:
%   G - graph with weights corresponding to the distances between each pair
%   of points if there is an edge between them

N = size(points,1);
W = sparse(N,N);
for i = 1:N
    for j = i+1:N
        %[dist,~] = distance(points(i,1),points(i,2),points(j,1),points(j,2),wgs84Ellipsoid);
        dist = norm(points(i,:)-points(j,:));
        distNoise = rand()*2*eta*dist - eta*dist; %uniform distribution in range [-eta*dist,eta*dist]
        dist = dist+distNoise;
        if(dist <= rho)
            W(i,j) = dist;
            W(j,i) = dist;
        end
    end
end
G = graph(W);
        

end


model = stlread('dazuoye.stl');
points = model.Points;
T = model.ConnectivityList;
trisurf(T,points(:,1),points(:,2),points(:,3),'FaceColor','none');
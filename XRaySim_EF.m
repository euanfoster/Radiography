clear all

%%
%Defining geometry

tic

%Defining the Steel Block Bounding Box
X1 = [107.3; 132.7; 132.7; 107.3; 107.3; 132.7; 132.7; 107.3];  
Y1 = [50; 50; -50; -50; 50; 50; -50; -50];
Z1 = [50; 50; 50; 50; -50; -50; -50; -50];

%Defining Hole in Component
x_cyl = 120;
y_cyl = 0;
r_cyl = 2.5;
th = 0:pi/50:2*pi;
X2 = r_cyl * cos(th) + x_cyl;
Y2 = r_cyl * sin(th) + y_cyl;
%Getting arrays in correct format for alphaShape Command
Z2 = [true(length(X2));zeros(length(X2))];
Z2 = Z2(:,1);
Z2 = Z2*50;
X2 = [X2,X2];
X2 = X2.';
Y2 = [Y2,Y2];
Y2 = Y2.';

%Creating alphaShape of component
plate = alphaShape(X1,Y1,Z1);
hole = alphaShape(X2,Y2,Z2);

%Plotting component
figure(01)
clf
plot(plate)
hold on
plot(hole)
axis equal
title('Steel Plate with Hole')
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z axis')
%%
%Dicretising Detector and Ray Paths

%Calculating FOV and detector size
r_block = sqrt(max(Y1)^2 + max(Z1)^2);
FOV = atan(r_block/(min(X1)));
r_detector = round(240*tan(FOV));

%Discretising the detector in Cartesian Coordinates
prompt = {'Enter Detector Pixels in Y Direction:','Enter Detector Pixels in Z Direction:'};
title = 'Input';
dims = [1 35];
definput = {'500','500'};
answer = inputdlg(prompt,title,dims,definput);
pix_y = str2double(answer(1,1)); %720;
pix_z = str2double(answer(2,1)); %1280;

detector_y = round(r_detector);
detector_z = round(r_detector*pix_z/pix_y);
detector_Y = linspace(-detector_y,detector_y,pix_y)';
detector_Z = linspace(-detector_z,detector_z,pix_z);
x = 240;

%Spherical Coordinates at the rectangulate detector
%projecting a rectangular beam through space along x axis
detector_R = round(sqrt(x.^2 + detector_Y.^2 + detector_Z.^2),0);
phi = atan(detector_Y./detector_Z);
theta = asin(detector_Z./(detector_R.*cos(phi)));

%creating 3D array of cartesian coordinates of the rectangular beam/array
inc = 0.1;
X_R = ones(pix_y,pix_z,240/inc);
for i = 1:size(X_R,3)
    X_R(:,:,i)= X_R(:,:,i)*inc*i;
end
Z_R = round((X_R./cos(theta)).*sin(theta).*cos(phi),0);
Y_R = round((X_R./cos(theta)).*sin(theta).*sin(phi)*-1,0);


%%
%Finding the coordinates of the rays that are in the block of steel
%Storing the ray's max/min coordinates in the block of steel and hole
%Calculating length of each ray in the steel and air


%Boolean 3D array of coordinates in and outside the block of steel
tf_plate = inShape(plate,X_R,Y_R,Z_R);
tf_hole = inShape(hole,X_R,Y_R,Z_R);

%Finding & storing the min indexes of the array when rays leave and enter the block of steel
plate_ind_min = zeros(pix_y,pix_z);
summation_plate = zeros(pix_y,pix_z);
for i = 1:pix_z
    for j = 1:pix_y
           summation_plate(j,i) = round(sum(tf_plate(j,i,:)),0);
           if summation_plate(j,i) ~= 0
               plate_ind_min(j,i) = find(tf_plate(j,i,:),1,'first');
           else
               plate_ind_min(j,i) = 0;
           end
    end
end

%Finding & storing the max indexes of the array when rays leave and enter the block of steel
plate_ind_max = zeros(pix_y,pix_z);
for i = 1:pix_z
    for j = 1:pix_y
           if summation_plate(j,i) ~= 0
               plate_ind_max(j,i) = find(tf_plate(j,i,:),1,'last');
           else
               plate_ind_max(j,i) = 0;
           end
    end
end

%Finding storing X, Y & Z Coordinates of rays entering and leaving the block of steel
X_Max = zeros(pix_y,pix_z);
Y_Max = zeros(pix_y,pix_z);
Z_Max = zeros(pix_y,pix_z);
X_Min = zeros(pix_y,pix_z);
Y_Min = zeros(pix_y,pix_z);
Z_Min = zeros(pix_y,pix_z);

%Finding max coorindates of ray in the block of steel
for i = 1:pix_z
    for j = 1:pix_y
            if plate_ind_max(j,i) ~= 0
                X_Max(j,i) = X_R(j,i,plate_ind_max(j,i));
                Y_Max(j,i) = Y_R(j,i,plate_ind_max(j,i));
                Z_Max(j,i) = Z_R(j,i,plate_ind_max(j,i));
            else
                X_Max(j,i) = 0;
                Y_Max(j,i) = 0;
                Z_Max(j,i) = 0;
            end
    end
end

%Finding min coordinates of ray in the block of steel
for i = 1:pix_z
    for j = 1:pix_y
            if plate_ind_min(j,i) ~= 0
                X_Min(j,i) = X_R(j,i,plate_ind_min(j,i));
                Y_Min(j,i) = Y_R(j,i,plate_ind_min(j,i));
                Z_Min(j,i) = Z_R(j,i,plate_ind_min(j,i));
            else
                X_Min(j,i) = 0;
                Y_Min(j,i) = 0;
                Z_Min(j,i) = 0;
            end
    end
end

%Calculating length of each ray in the block of steel
deltax = X_Max-X_Min;
deltay = Y_Max-Y_Min;
deltaz = Z_Max-Z_Min;
length_steel = sqrt(deltax.^2 + deltay.^2 + deltaz.^2);

%Clearing variables for the hole
clear('X_Max','Y_Max','Z_Max','X_Min','Y_Min','Z_Min','deltax','deltay','deltaz')

%Finding & storing the min indexes of the array when rays leave and enter the hole
hole_ind_min = zeros(pix_y,pix_z);
summation_hole = zeros(pix_y,pix_z);
for i = 1:pix_z
    for j = 1:pix_y
           summation_hole(j,i) = round(sum(tf_hole(j,i,:)),0);
           if summation_hole(j,i) ~= 0
               hole_ind_min(j,i) = find(tf_hole(j,i,:),1,'first');
           else
               hole_ind_min(j,i) = 0;
           end
    end
end

%Finding & storing the max indexes of the array when rays leave and enter the hole
hole_ind_max = zeros(pix_y,pix_z);
for i = 1:pix_z
    for j = 1:pix_y
           if summation_hole(j,i) ~= 0
               hole_ind_max(j,i) = find(tf_hole(j,i,:),1,'last');
           else
               hole_ind_max(j,i) = 0;
           end
    end
end

%Finding storing X, Y & Z Coordinates of rays entering and leaving the hole
X_Max = zeros(pix_y,pix_z);
Y_Max = zeros(pix_y,pix_z);
Z_Max = zeros(pix_y,pix_z);
X_Min = zeros(pix_y,pix_z);
Y_Min = zeros(pix_y,pix_z);
Z_Min = zeros(pix_y,pix_z);

%Finding max coorindates of ray in the hole
for i = 1:pix_z
    for j = 1:pix_y
            if hole_ind_max(j,i) ~= 0
                X_Max(j,i) = X_R(j,i,hole_ind_max(j,i));
                Y_Max(j,i) = Y_R(j,i,hole_ind_max(j,i));
                Z_Max(j,i) = Z_R(j,i,hole_ind_max(j,i));
            else
                X_Max(j,i) = 0;
                Y_Max(j,i) = 0;
                Z_Max(j,i) = 0;
            end
    end
end

%Finding min coorindates of ray in the hole
for i = 1:pix_z
    for j = 1:pix_y
            if hole_ind_min(j,i) ~= 0
                X_Min(j,i) = X_R(j,i,hole_ind_min(j,i));
                Y_Min(j,i) = Y_R(j,i,hole_ind_min(j,i));
                Z_Min(j,i) = Z_R(j,i,hole_ind_min(j,i));
            else
                X_Min(j,i) = 0;
                Y_Min(j,i) = 0;
                Z_Min(j,i) = 0;
            end
    end
end

%Calculating length of each ray in the hole
deltax = X_Max-X_Min;
deltay = Y_Max-Y_Min;
deltaz = Z_Max-Z_Min;
length_hole = sqrt(deltax.^2 + deltay.^2 + deltaz.^2);

%clearing variables
clear('X_Max','Y_Max','Z_Max','X_Min','Y_Min','Z_Min','deltax','deltay','deltaz')

%calculating the length of each ray if it passes through the hole
length_steel = length_steel - length_hole;
%calculating the length of each reay in air
length_air = detector_R - length_steel;

%%
%Calculating the total attentuation and plotting result
%attenuation coefficients
alpha_air = 6.358e-04.*0.001225;
alpha_steel = 7.87.*5.995e-04;

%calculating the attenuation of the rays in air and steel
att_air = exp(-1.*alpha_air.*length_air);
att_steel = exp(-1.*alpha_steel.*length_steel);

%calculating total image
att_total = (att_air + att_steel);

toc

figure(02)
clf
pcolor(Z_R(:,:,2400),Y_R(:,:,2400),att_total)
set(gca, 'clim', [1.85,2]);
shading interp;
colormap gray
MOP = colormap;
MOP = flipud(MOP);
colormap(MOP);
colorbar


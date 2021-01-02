%% script to read abaqus report and plot slip traces
%Yi Xiong
%University of Oxford
%Sep 2020

% Need to convert .rpt into .xlsx using space as delimiter
clear all;
fid='bottom.xlsx';    %input file name 
f=xlsread(fid);
element_num=f(:,1);
if_element=element_num>1;        
index=[];
% find rows which contain gnd data
for i=1:length(if_element)
    if if_element(i)==1;
        index=[index;i];
    end
end
data=[];
% store the useful data
for i=1:length(index)
    data(i,:)=f(index(i),:);
end

% Get slip direction and plane normals for 3 prism slip systems
[slip_dir,slip_dnor]=slip_systems;
prismdir=slip_dir(4:6,:);
prismdnor=slip_dnor(4:6,:);

% Define grain orientations, extract SDV1-SDV9 to form rotation matrix
for i=1:length(data(:,1))
    R(1,1,i)=data(i,2);
    R(1,2,i)=data(i,3);
    R(1,3,i)=data(i,4);
    R(2,1,i)=data(i,5);
    R(2,2,i)=data(i,6);
    R(2,3,i)=data(i,7);
    R(3,1,i)=data(i,11);
    R(3,2,i)=data(i,12);
    R(3,3,i)=data(i,13);
end
% Find the slip system with highest gnd
trace(:,1:3)=abs(data(:,8:10));        %read gnd values for 3 prism slip 
e=[0,0,1];                             %sample surface plane normal
for i=1:length(trace(:,1))
    gnd=abs(trace(i,1:3));
    M(i)=max(gnd);                     %Find the active slip system
    if M(i)==0
        I(i)=1;
    else
        I(i)=find(gnd==M(i));          %Get index of the active slip system
    end
    trace(i,4:6)=prismdnor(I(i),:);      %slip normal vector of the slip system with highest gnd
    norm(i,1:3)=R(:,:,i)*trace(i,4:6)';  %rotate into sample coordinates      
    line(i,1:3)=cross(e,norm(i,1:3));    %cross product of surface normal and the slip plane normal
    line(i,2)=-line(i,2);                %flip along x-axis, so that to match the EBSD coordinate and MATLAB coordinate
end

column=139;           %number of element in one column of the model

% calculate the coordinate of each element
for i=1:length(data(:,1))
    yvalues(i,1)=fix(data(i,1)/column);
    xvalues(i,1)=mod(data(i,1),column);
    if xvalues(i,1)==0
        xvalues(i,1)=column;
    end
end
yval=column-xvalues+1;
xval=yvalues-(min(yvalues)-1);

for i=1:length(xval)
    l=M(i)/mean(M(i));       %define the length of the line segments to be proportional to gnd values   
    plot([xval(i)-l*line(i,1) xval(i)+l*line(i,1)],[yval(i)-l*line(i,2) yval(i)+l*line(i,2)],'r');   %plot the direction of the lines
    hold on 
end






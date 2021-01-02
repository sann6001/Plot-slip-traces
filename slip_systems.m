function [dir,dnor] = slip_systems()
%% Slipsystems
a = 1; c =1.587; % 1.6236; % scale 3 component of cartiesian vector by  c/a for specific material 
a1 = a*[-1/2 sqrt(3)/2 0];
a2 = a*[-1/2 -sqrt(3)/2 0];
a3 = a*[1 0 0];
a4 = c*[0 0 1];

% 3<a> <11-20> slip directions each on 4 planes (Ba, Pr, PyI, PyII)
s(1,:) = 1/3*[2 -1 -1 0]; %b1
s(2,:) = 1/3*[-1 -1 2 0]; %b2
s(3,:) = 1/3*[-1 2 -1 0]; %b3
%6 <c+a> <11-23>  slip directions each on 4 planes (Pr, PyI, PyII, sP)
s(4,:)=1/3*[2 -1 -1 3]; % b4
s(5,:)=1/3*[2 -1 -1 -3]; % b5
s(6,:)=1/3*[-1 -1 2 3]; % b6
s(7,:)=1/3*[-1 -1 2 -3]; % b7
s(8,:)=1/3*[-1 2 -1 3]; % b8
s(9,:)=1/3*[-1 2 -1 -3]; % b9 


for i =1:9
    S(i,:) =s(i,1)*a1 + s(i,2)*a2 + s(i,3)*a3+s(i,4)*a4;       
  
    %S(i,:) = S(i,:)/norm(S(i,:));
end

% slip plane normal for Basal plane (1Ba, 2Ba, 3Ba)
n = [0 0 0 1];
Ba = a4/n(4);

% slip plane normals for 6 1st order pyramidal {10-11} (PyI,PyII)
n = [0 1 -1 1]; % 1PyI = 6PyII = 9PyI
PyI(1,:) = cross(a2/n(2) - a3/n(3), a4/n(4) - a3/n(3));

n = [1 -1 0 1]; % 2PyI = 5PyII = 8PyI
PyI(2,:) = cross(a1/n(1) - a2/n(2), a4/n(4) - a2/n(2));

n = [-1 0 1 1]; % 3PyI = 4PyI = 7PyI
PyI(3,:) = cross(a4/n(4) - a3/n(3), a1/n(1) - a3/n(3));

n = [0 -1 1 1]; % 1PyII = 7PyII = 8PyII
PyII(1,:) = cross(a2/n(2) - a3/n(3), a4/n(4) - a3/n(3));

n = [-1 1 0 1]; % 2PyII = 4PyII = 9PyII
PyII(2,:) = cross(a1/n(1) - a2/n(2), a4/n(4) - a2/n(2));

n= [1  0 -1 1]; % 3PyII = 5Py1 = 6PyI
PyII(3,:) = cross(a4/n(4) - a3/n(3), a1/n(1) - a3/n(3));

PyI(4,:) = PyI(3,:);
PyI(5,:) = PyII(3,:);
PyI(6,:) = PyII(3,:);
PyI(7,:) = PyI(3,:);
PyI(8,:) = PyI(2,:);
PyI(9,:) = PyI(1,:);

% slip plane normals for 6 2nd order pyramidal {11-22} (sP)
n = [-2 1 1 2]; % 4sP
sP(4,:) = cross(a3/n(3) - a2/n(2), a4/n(4) - a2/n(2));

n = [2 -1 -1 2]; % 5sP
sP(5,:) = cross(a3/n(3) - a2/n(2), a4/n(4) - a2/n(2));

n = [1 1 -2 2]; % 6sP
sP(6,:) = cross(a4/n(4) - a2/n(2), a1/n(1) - a2/n(2));

n = [-1 -1 2 2]; % 7sP
sP(7,:) = cross(a4/n(4) - a2/n(2), a1/n(1) - a2/n(2));


n = [1 -2 1 2]; % 8sP
sP(8,:) = cross(a4/n(4) - a1/n(1), a3/n(3) - a1/n(1));

n = [-1 2 -1 2]; % 9sP
sP(9,:) = cross(a4/n(4) - a1/n(1), a3/n(3) - a1/n(1));


% slip plane normals for 6 prismatic {10-10} planes (Pr)
 n = [0 1 -1 0]; % 1Pr = 4Pr = 5Pr
 Pr(1,:) = cross(a2/n(2)-a3/n(3),a4);
  
 n = [1 -1 0 0]; % 2Pr = 6Pr = 7Pr
 Pr(2,:) = cross(a1/n(1)-a2/n(2),a4);

 n = [-1 0 1 0]; % 3Pr = 8Pr = 9Pr note type in 9Pr normal
 Pr(3,:) = cross(a3/n(3)-a1/n(1),a4 ) ;

% plane normals or {12-31} planes (P121)

 n=[1 2 -3 1];
 P121(1,:)= cross(a1/n(1)-a3/n(3), a4/n(4)-a3/n(3));
 
 n=[-1 -2 3 1];
 P121(2,:)= cross(a1/n(1)-a3/n(3), a4/n(4)-a3/n(3));
 
 n=[2 1 -3 1];
 P121(3,:)= cross(a1/n(1)-a3/n(3), a4/n(4)-a3/n(3));
 
 n=[-2 -1 3 1];
 P121(4,:)= cross(a1/n(1)-a3/n(3), a4/n(4)-a3/n(3));
  
% normalise vectors
for i=1:3;
    dir(i,:)=S(i,:)/norm(S(i,:));
    dir(i+3,:)=S(i,:)/norm(S(i,:));
end
for i=4:9;
    dir(i+3,:)=S(i,:)/norm(S(i,:));
end
for i=1:3;
    dnor(i,:)=Ba/norm(Ba);
end
for i=4:6;
    dnor(i,:)=Pr(i-3,:)/norm(Pr(i-3,:));
end
for i=7:12;
    dnor(i,:)=PyI(i-3,:)/norm(PyI(i-3,:));
end
for i=13:18;
    dnor(i,:)=sP(i-9,:)/norm(sP(i-9,:));
end
 N = sP(4,:);
 A = a1/n(1);
 B = a2/n(2);
 C = a3/n(3);
 
 for i=1:4;
     N121(i,:)= P121(i,:)/norm(P121(i,:));
 end

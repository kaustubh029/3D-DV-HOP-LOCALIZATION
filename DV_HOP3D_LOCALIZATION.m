%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3D DV-Hop Localization~~~~~~~~~~~~~~~~~~~~~~~~
clear,close all;
clc;

maxTrial = 5;
count = 0;
BorderLength = 50;                                         % Length of area of the network
NodeAmount = 250;                                          % Total number of sensor nodes
BeaconAmount = 25;                                         % Number of Anchor nodes
UNAmount = NodeAmount - BeaconAmount;                      % Number of Unknown sensor nodes
Range = 15;                                                % Communication range of sensor nodes
RangeError = .2;                                           % Communication ranging error


DistanceAll = zeros(NodeAmount,NodeAmount);                % Distance between node to node
hops = zeros(NodeAmount,NodeAmount);
hopsize = zeros(BeaconAmount,1);
DistanceAtUN = zeros(BeaconAmount,UNAmount);               % Distance between anchor nodes to unknown nodes
A = zeros(BeaconAmount-1,3);                               % A matrix of least mean square method
B = zeros(BeaconAmount-1,1);                               % B matrix of least mean square method
X = zeros(3,UNAmount);                                     % X matrix of least mean square method
coordError = zeros(1,UNAmount);                            % Difference(distance) between estimated and actual coordinates of unknown nodes
error = zeros(maxTrial,1);                                 % Error in calculation of unknown nodes


for trial = 1:maxTrial
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Coordinates = BorderLength.*rand(3,NodeAmount);                                                                                                         % Coordinates of nodes
    RandomRangeError = rand(NodeAmount,NodeAmount) * Range * RangeError;
    BeaconCoords = [Coordinates(1,1:BeaconAmount); Coordinates(2,1:BeaconAmount); Coordinates(3,1:BeaconAmount)];                                           % Coordinates of anchor nodes
    UNCoords = [Coordinates(1,(BeaconAmount+1):NodeAmount); Coordinates(2,(BeaconAmount+1):NodeAmount); Coordinates(3,(BeaconAmount+1):NodeAmount)];        % Coordinates of unknown nodes
    
    plot3(Coordinates(1,1:BeaconAmount), Coordinates(2,1:BeaconAmount), Coordinates(3,1:BeaconAmount), 'r*', Coordinates(1,(BeaconAmount+1):NodeAmount), Coordinates(2,(BeaconAmount+1):NodeAmount), Coordinates(3,(BeaconAmount+1):NodeAmount), 'k.')
    xlim([0,BorderLength]);
    ylim([0,BorderLength]);
    zlim([0,BorderLength]);
    title('Random Distribution of Anchor and Sensor nodes')

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i = 1:NodeAmount
        for j = 1:NodeAmount
            DistanceAll(i,j) = ((Coordinates(1,i)-Coordinates(1,j))^2 + (Coordinates(2,i)-Coordinates(2,j))^2 + (Coordinates(3,i)-Coordinates(3,j))^2)^0.5;
            if (DistanceAll(i,j) <= Range - RandomRangeError(i,j)) && (DistanceAll(i,j) > 0)
                hops(i,j)=1;
            elseif i==j
                hops(i,j)=0;
            else
                hops(i,j)=inf;
            end
         end
    end

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for k = 1:NodeAmount
        for i = 1:NodeAmount
            for j = 1:NodeAmount

                if hops(i,k) + hops(k,j) < hops(i,j)
                    hops(i,j) = hops(i,k) + hops(k,j);

                end
             end
        end
    end

    %~~~~~~~~~~~~~~~~~~~~~~~~~Hop-size by DV_HOP Method~~~
    hopsAtA = hops(1:BeaconAmount, 1:BeaconAmount);                 % Minimum number of hops between Anchor node to Anchor node
    DistanceAtA = DistanceAll(1:BeaconAmount, 1:BeaconAmount);      % Distance between Anchor node to Anchor node
    
    %tic

    for i = 1:BeaconAmount
        hopsize(i,1) = sum(DistanceAtA(i,:)) / sum(hopsAtA(i,:));
    end

    hopsAtUN = hops(1:BeaconAmount,(BeaconAmount+1):NodeAmount);    % Minimum number of hops between Anchor node to Unknown node

    %~~~~~~~~~~~~~~~~~~~~~~~~~~Distance estimation between Unknown node to Anchor node by DV_Hop method~~~
    for i = 1:UNAmount
        for j = 1:BeaconAmount
            [minhop,mini] = min(hopsAtUN(:,i));
            hopnum = hopsAtUN(j,i);
            hop = hopsize(mini,1);
            DistanceAtUN(j,i) = hopsize(j,1) * hopnum;
        end
    end

    %~~~~~~~~~~~~~~~~~~~~~~~~~ Matrix "A" of DV_Hop algorithm ~~~
    for i = 1:3
        for j = 1:(BeaconAmount-1)
            A(j,i) = 2 * (BeaconCoords(i,j)-BeaconCoords(i,BeaconAmount));
        end
    end

    %~~~~~~~~~~~~~~~~~~~~~~~~~~Matrix "B" of DV-Hop algorithm and Coordinate estimation of Unknown node~~~
    for m = 1:UNAmount 
        for i = 1:(BeaconAmount-1)
            B(i,1) = BeaconCoords(1,i)^2-BeaconCoords(1,BeaconAmount)^2+BeaconCoords(2,i)^2-BeaconCoords(2,BeaconAmount)^2+BeaconCoords(3,i)^2-BeaconCoords(3,BeaconAmount)^2+DistanceAtUN(BeaconAmount,m)^2-DistanceAtUN(i,m)^2;
        end
               X1 = inv(A'*A)*A'*B;
               X(1,m) = X1(1,1);      % estimated x coordinate of unknown node
               X(2,m) = X1(2,1);      % estimated y coordinate of unknown node
               X(3,m) = X1(3,1);      % estimated z coordinate of unknown node

     end
     %t = toc

     figure;
     hold on
     plot3(Coordinates(1,1:BeaconAmount), Coordinates(2,1:BeaconAmount), Coordinates(3,1:BeaconAmount),'ro')
     plot3(Coordinates(1,(BeaconAmount+1):NodeAmount), Coordinates(2,(BeaconAmount+1):NodeAmount), Coordinates(3,(BeaconAmount+1):NodeAmount),'k*', X(1,:), X(2,:), X(3,:), 'm.')
     for i = (BeaconAmount+1):NodeAmount
         plot3([Coordinates(1,i),X(1,i-BeaconAmount)],[Coordinates(2,i),X(2,i-BeaconAmount)],[Coordinates(3,i),X(3,i-BeaconAmount)],'b-')
     end
      title ('Anchor nodes (o), Sensor nodes (*) and estimated Position of sensor nodes (.)')
    hold off

    % ~~~~~~~~~~~~~~~~~~~~~~~~~Error calculation~~~~~~~~~~~~~~~~~~~~~~~~~
    for i = 1:UNAmount
        coordError(1,i) = (((X(1,i) - UNCoords(1,i))^2 + (X(2,i) - UNCoords(2,i))^2 + (X(3,i) - UNCoords(3,i))^2)^0.5);   
    end

    error(trial,1) = sum(coordError(1,:))/(UNAmount * Range);

    if isnan(error(trial,1))
        error(trail,1) = 0;
        count = count+1;
    end

    figure;
    plot (sort(coordError),'r*')
    title ('Error')

end
avgError = sum(error(:,1)) / (maxTrial-count)


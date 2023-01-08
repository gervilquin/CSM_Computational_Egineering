function Be = body_forces(X,g)
    
    Be = zeros(length(X(:,1)),3);

    for i =1:length(X(:,1))
        Be(i,:) = [g i 3];
    end
    
end
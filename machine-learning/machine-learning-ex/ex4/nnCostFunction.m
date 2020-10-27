function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%

% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));

% ====================== YOUR CODE HERE ======================
% Instructions: You should complete the code by working through the
%               following parts.
%
% Part 1: Feedforward the neural network and return the cost in the
%         variable J. After implementing Part 1, you can verify that your
%         cost function computation is correct by verifying the cost
%         computed in ex4.m
%
% Part 2: Implement the backpropagation algorithm to compute the gradients
%         Theta1_grad and Theta2_grad. You should return the partial derivatives of
%         the cost function with respect to Theta1 and Theta2 in Theta1_grad and
%         Theta2_grad, respectively. After implementing Part 2, you can check
%         that your implementation is correct by running checkNNGradients
%
%         Note: The vector y passed into the function is a vector of labels
%               containing values from 1..K. You need to map this vector into a 
%               binary vector of 1's and 0's to be used with the neural network
%               cost function.
%
%         Hint: We recommend implementing backpropagation using a for-loop
%               over the training examples if you are implementing it for the 
%               first time.
%
% Part 3: Implement regularization with the cost function and gradients.
%
%         Hint: You can implement this around the code for
%               backpropagation. That is, you can compute the gradients for
%               the regularization separately and then add them to Theta1_grad
%               and Theta2_grad from Part 2.
%

ohey = [y == 1:10]; % one-hot encoding, return y matrix 5000*10

Y = zeros(m,num_labels);
for i = 1:m
    Y(i, y(i)) = 1;
end

% Part #1
% solution #1 (for loop)
for i = 1:m
    a1 = sigmoid([1,X(i,:)]*Theta1');
    a2 = sigmoid([1,a1]*Theta2');
    for j = 1:num_labels
        if j == y(i)
            J = J -log(a2(j)) ; %notice that y is a matrix only contains 0 or 1
        else
            J = J - log(1-a2(j));
        end
    end
end
J = J/m;

% % solution #2 (regularized version)
% hx = sigmoid([ones(m,1),sigmoid([ones(m,1),X] * Theta1')] * Theta2');
% J = 1/m * sum(sum(-ohey.*log(hx)-(1-ohey).*log(1-hx)));

% Regularized
%J = J + 0.5*lambda/m*( sum(sum(Theta1(:, 2:input_layer_size+1).^2))+  sum(sum(Theta2(:,2:hidden_layer_size+1).^2)));
J = J + 0.5*lambda./m.*( sum(sum(Theta1(:, 2:end).^2))+  sum(sum(Theta2(:,2:end).^2)));

% Backpropagration
for t=1:m
    a_1 = X(t,:); % 1*400, first
    z_2 = [1,a_1]*Theta1'; % 1*25
    a_2 = sigmoid(z_2); %1*25, second
    z_3 = [1,a_2]*Theta2'; %1*10
    a_3 = sigmoid(z_3); %1*10, third output
    
    %delta_3 = a_3 - ohey(t,:); %1*10 !!!!!!!!!!!!!!!!!!!every time I use
    %logical array, there is wierd error 'Matrix dimension'
    
    delta_3 = a_3 - Y(t,:); %1*10
    delta_2 = (delta_3*Theta2(:,2:end)).*sigmoidGradient(z_2); %25
    Theta2_grad = Theta2_grad + delta_3'*[1,a_2]; %10*25
    Theta1_grad = Theta1_grad + delta_2'*[1,a_1]; %25*401
end

Theta2_grad = Theta2_grad / m; 
Theta1_grad = Theta1_grad / m;

Theta1_grad(:,2:end) = Theta1_grad(:,2:end) + lambda / m * Theta1(:, 2:end);
Theta2_grad(:,2:end) = Theta2_grad(:,2:end) + lambda / m * Theta2(:, 2:end);

% Unroll gradients

grad = [Theta1_grad(:) ; Theta2_grad(:)];





% -------------------------------------------------------------

% =========================================================================

% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:)];


end

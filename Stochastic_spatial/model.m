%% 2D lattice of states of the system
%% a grid of boolean states 

SCS = rand(200);
disp(SCS)
prob_birth_coeff = 0.55; %% wiggle later from 0 to 1


T=100;
%initialise


for i=1:T
    SCS=updateSCS(SCS, prob_birth_coeff);
    %if mod(i, 2) == 0 
        imshow(SCS);
    %end
    disp(size(SCS));
end



function ret=updateSCS(SCS, prob_birth_coeff) %for each point find number of neighbours with state FALSE, if FASLE<

    matrix_mask = [0.707, 1, 0.707; 1, 0, 1; 0.707, 1, 0.707];
    neighbours_matrix = conv2(SCS, matrix_mask);
    neighbours_matrix = neighbours_matrix(2:length(SCS)+1, 2:length(SCS));
    neighbours_proportion = neighbousrs_matrix / 8; %% change neighbours to proportions for probabilities
 
    birth_prob_mat = neighbours_proportion;
    noisy_norm_prob_mat = birth_prob_mat + (0.03*unifrnd(200,199));

    ret = noisy_norm_prob_mat < prob_birth_coeff;

end
